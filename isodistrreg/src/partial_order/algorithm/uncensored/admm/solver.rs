//! The ADMM iteration for `min 0.5 u^T u + q^T u  s.t.  A u >= 0`.
//!
//! This is the operator-splitting scheme of Stellato et al. (2020) specialized to the
//! partial-order quadratic program. Two properties of that program do the work:
//!
//!  * `P = I` in energy coordinates (see the caller), so eliminating the slack from the
//!    indefinite ADMM system leaves the normal equations `K = (1 + sigma) I + rho A^T A`,
//!    an SPD matrix with `lambda_min >= 1 + sigma`. Plain Cholesky suffices, and the
//!    condition number is known in closed form.
//!  * The pattern of `K` is the cover graph plus a diagonal and never changes across the
//!    threshold sequence, so the symbolic factorization is computed once per fit and only
//!    values are refreshed -- on a `rho` change, or on the hazard-rate path where `A`'s
//!    values are rewritten per threshold.
//!
//! The cone is the nonnegative orthant, so the projection is `max(0, .)` and the entries
//! it clips are, bit for bit, the active constraint set. `Workspace::clipped` records them
//! for the polishing step.
use super::constraints::Constraints;
use super::kkt::KktFactor;

/// Regularization added to `P` in the ADMM system. OSQP's default.
pub(crate) const SIGMA: f64 = 1e-6;
/// Over-relaxation factor. OSQP's default.
pub(crate) const ALPHA: f64 = 1.6;
/// Initial penalty. OSQP's default.
pub(crate) const RHO_INITIAL: f64 = 0.1;
/// Penalty bounds, matching `OSQP_RHO_MIN` / `OSQP_RHO_MAX`.
pub(crate) const RHO_MIN: f64 = 1e-6;
pub(crate) const RHO_MAX: f64 = 1e6;
/// A new `rho` is adopted only if it differs by more than this factor, because adopting
/// one costs a numeric refactorization. Matches `OSQP_ADAPTIVE_RHO_TOLERANCE`.
pub(crate) const ADAPTIVE_RHO_TOLERANCE: f64 = 5.0;
/// Iterations between residual checks, and hence between `rho` updates.
///
/// A check costs two matrix-vector products, `O(m)`, against an iteration's triangular
/// solve at `O(nnz(L))`; on these problems that is a tenth to a quarter of an iteration.
/// Checking more often therefore pays for itself by not overshooting: consecutive
/// thresholds are warm-started and often converge in far fewer iterations than the
/// interval, and every iteration past convergence is pure waste. Measured on the real
/// instances, moving from 25 to 10 cut iterations by 9-18%; going below 10 gained nothing.
pub(crate) const CHECK_INTERVAL: u32 = 10;

/// Relative forward error of an f32 solve of `K x = b`, per unit of `cond(K)`.
///
/// Measured on this system across `cond(K)` from 1e2 to 1e6: the *backward* error is flat
/// at 0.26-0.58 f32 ULP, so the factorization itself is backward stable and every bit of
/// forward error is conditioning. The observed ratio is 7.2e-9, 6.1e-9, 7.2e-9 at
/// `cond` 1e2, 1e4, 1e6 -- about 0.06 f32 ULP per unit of condition number.
pub(crate) const F32_ERROR_PER_COND: f64 = 7.2e-9;

/// Relative forward error of an f64 solve, per unit of `cond(K)`.
///
/// The classical bound for a backward-stable solve, and what the measured f64 residual
/// floor tracks: at `cond(K) = 3.6e6` the dual residual cannot be pushed below about
/// 8e-10, which is this constant times that condition number.
pub(crate) const F64_ERROR_PER_COND: f64 = f64::EPSILON;

/// Fraction of `eps_abs` the inexact `x`-update is allowed to consume.
///
/// ADMM tolerates an inexact linear solve, but not one comparable to the tolerance it is
/// trying to reach: the residuals would plateau above `eps_abs` and every solve would burn
/// its whole iteration budget without converging. Keeping the solve error an order of
/// magnitude below the target leaves the termination test meaningful.
pub(crate) const F32_ERROR_MARGIN: f64 = 0.1;

/// Why the iteration stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Status {
    /// Both residuals met their tolerance.
    Solved,
    /// The iteration budget ran out. `Workspace::z` is still feasible by construction.
    MaxIterations,
}

/// Tuning handed to a single solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Settings {
    pub(crate) verbose: bool,
    pub(crate) eps_abs: f64,
    pub(crate) eps_rel: f64,
    pub(crate) max_iter: u32,
}

/// What one solve did, for diagnostics and for the caller's convergence accounting.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Outcome {
    pub(crate) status: Status,
    pub(crate) iterations: u32,
    pub(crate) refactorizations: u32,
    /// Primal residual `||A u - z||_inf` at the returned iterate.
    pub(crate) primal_residual: f64,
    /// Dual residual `||u + q + A^T y||_inf` at the returned iterate.
    pub(crate) dual_residual: f64,
    /// The penalty the solve ended on, carried into the next threshold.
    pub(crate) rho: f64,
    /// Whether the final factorization was stored in f32.
    pub(crate) used_f32: bool,
}

/// Iterate state plus every scratch buffer, reused across the threshold sequence.
///
/// `x`, `z` and `y` persist between solves: consecutive thresholds differ only in `q`,
/// so carrying the iterate forward is the warm start.
pub(crate) struct Workspace {
    pub(crate) x: Vec<f64>,
    pub(crate) z: Vec<f64>,
    pub(crate) y: Vec<f64>,
    /// `clipped[r]` iff the projection zeroed row `r` on the final update, i.e. the
    /// active set. Recorded rather than recovered from `z` so that a row which happens
    /// to land on exactly zero without being clipped is not misread as active.
    pub(crate) clipped: Vec<bool>,
    /// `A x` at the returned iterate; the polish needs it and it is already computed.
    pub(crate) au: Vec<f64>,
    rhs: Vec<f64>,
    aty: Vec<f64>,
    relaxed: Vec<f64>,
    column_scratch: Vec<f64>,
    rho: f64,
    /// `rho` the current factor was built with; `f64::NAN` when there is no factor yet.
    factor_rho: f64,
    /// Precision the current factor was built in.
    factor_used_f32: bool,
}

impl Workspace {
    pub(crate) fn new(n: usize, m: usize) -> Self {
        Self {
            x: vec![0.0; n],
            z: vec![0.0; m],
            y: vec![0.0; m],
            clipped: vec![false; m],
            au: vec![0.0; m],
            rhs: vec![0.0; n],
            aty: vec![0.0; n],
            relaxed: vec![0.0; m],
            column_scratch: vec![0.0; n],
            rho: RHO_INITIAL,
            factor_rho: f64::NAN,
            factor_used_f32: false,
        }
    }

    /// Force a refactorization on the next solve, because `A`'s values changed.
    pub(crate) fn invalidate_factor(&mut self) {
        self.factor_rho = f64::NAN;
    }
}

/// Infinity norm that *propagates* NaN.
///
/// `f64::max` returns the non-NaN operand, so the obvious fold silently drops a NaN and a
/// diverged iterate would report a zero residual and be accepted as solved. Propagating
/// instead makes every subsequent comparison false, so such an iterate can only ever be
/// reported as not converged.
fn norm_inf(v: &[f64]) -> f64 {
    let mut worst = 0.0f64;
    for &value in v {
        let magnitude = value.abs();
        if magnitude.is_nan() {
            return f64::NAN;
        }
        if magnitude > worst {
            worst = magnitude;
        }
    }
    worst
}

/// Solve `min 0.5 u^T u + q^T u  s.t.  A u >= 0`, leaving the primal iterate in `ws.x`.
///
/// `ws` carries the warm start in and the solution out. `factor` must have been built
/// with `KktFactor::new` for the same `n` and edge list; this routine refreshes its
/// numeric values whenever `rho` moves or `Workspace::invalidate_factor` was called.
pub(crate) fn solve(
    a: &Constraints<'_>,
    q: &[f64],
    settings: &Settings,
    ws: &mut Workspace,
    factor: &mut KktFactor,
) -> Outcome {
    debug_assert_eq!(q.len(), a.n);
    let lambda_max_bound = a.spectral_norm_squared_bound(&mut ws.column_scratch);
    let norm_q = norm_inf(q);
    let rho_limit = rho_ceiling(lambda_max_bound, settings.eps_abs);
    ws.rho = ws.rho.min(rho_limit);

    let mut refactorizations = 0u32;
    // `sqrt(eps)` is the smallest step whose square is still resolvable; below it the
    // residual ratio is noise and re-tuning `rho` would chase rounding.
    let residual_floor = f64::EPSILON.sqrt();

    refresh_factor(
        a,
        lambda_max_bound,
        settings.eps_abs,
        ws,
        factor,
        &mut refactorizations,
    );

    let mut status = Status::MaxIterations;
    let mut iterations = settings.max_iter;
    let mut primal_residual = f64::INFINITY;
    let mut dual_residual = f64::INFINITY;

    for iteration in 1..=settings.max_iter {
        // rhs = sigma x - q + A^T (rho z - y); reuse `w` to hold `rho z - y`.
        for ((scaled_r, &z_r), &y_r) in ws.relaxed.iter_mut().zip(&ws.z).zip(&ws.y) {
            *scaled_r = ws.rho * z_r - y_r;
        }
        a.mul_transpose(&ws.relaxed, &mut ws.aty);
        for ((rhs_c, &x_c), (&q_c, &aty_c)) in
            ws.rhs.iter_mut().zip(&ws.x).zip(q.iter().zip(&ws.aty))
        {
            *rhs_c = SIGMA * x_c - q_c + aty_c;
        }
        factor.solve_in_place(&mut ws.rhs);

        // `x_tilde` is now in `ws.rhs`, and the second block row of the ADMM system makes
        // `z_tilde` exactly `A x_tilde`, so no separate slack update is needed.
        a.mul(&ws.rhs, &mut ws.au);

        // x <- alpha x_tilde + (1 - alpha) x
        for (x_c, &xt_c) in ws.x.iter_mut().zip(&ws.rhs) {
            *x_c = ALPHA * xt_c + (1.0 - ALPHA) * *x_c;
        }
        // `relaxed` is the over-relaxed row value v = alpha z_tilde + (1 - alpha) z. The
        // projection acts on v + y/rho, but the dual step is driven by v itself: the
        // scaled-dual form of ADMM is u <- u + v - z, i.e. y <- y + rho (v - z). Feeding
        // the shifted quantity into the dual step instead would add y a second time and
        // the iteration would double every sweep.
        for ((relaxed_r, &zt_r), &z_r) in ws.relaxed.iter_mut().zip(&ws.au).zip(&ws.z) {
            *relaxed_r = ALPHA * zt_r + (1.0 - ALPHA) * z_r;
        }
        for (((z_r, clipped_r), &relaxed_r), &y_r) in
            ws.z.iter_mut()
                .zip(&mut ws.clipped)
                .zip(&ws.relaxed)
                .zip(&ws.y)
        {
            let shifted = relaxed_r + y_r / ws.rho;
            *clipped_r = shifted <= 0.0;
            *z_r = if shifted > 0.0 { shifted } else { 0.0 };
        }
        for ((y_r, &relaxed_r), &z_r) in ws.y.iter_mut().zip(&ws.relaxed).zip(&ws.z) {
            *y_r += ws.rho * (relaxed_r - z_r);
        }

        let last = iteration == settings.max_iter;
        if iteration % CHECK_INTERVAL != 0 && !last {
            continue;
        }

        // Residuals at the current iterate. `A x` is recomputed because `ws.au` holds
        // `A x_tilde`, which the relaxation step has since moved away from.
        a.mul(&ws.x, &mut ws.au);
        a.mul_transpose(&ws.y, &mut ws.aty);
        let mut r_prim = 0.0f64;
        for (&au_r, &z_r) in ws.au.iter().zip(&ws.z) {
            let residual = (au_r - z_r).abs();
            if residual.is_nan() {
                r_prim = f64::NAN;
                break;
            }
            if residual > r_prim {
                r_prim = residual;
            }
        }
        let mut r_dual = 0.0f64;
        for ((&x_c, &q_c), &aty_c) in ws.x.iter().zip(q).zip(&ws.aty) {
            let residual = (x_c + q_c + aty_c).abs();
            if residual.is_nan() {
                r_dual = f64::NAN;
                break;
            }
            if residual > r_dual {
                r_dual = residual;
            }
        }
        primal_residual = r_prim;
        dual_residual = r_dual;

        let scale_prim = norm_inf(&ws.au).max(norm_inf(&ws.z));
        let scale_dual = norm_inf(&ws.x).max(norm_inf(&ws.aty)).max(norm_q);
        let eps_prim = settings.eps_abs + settings.eps_rel * scale_prim;
        let eps_dual = settings.eps_abs + settings.eps_rel * scale_dual;

        if settings.verbose {
            eprintln!(
                "  admm iter {iteration:>5}  r_prim {r_prim:.3e} (<= {eps_prim:.3e})  \
                 r_dual {r_dual:.3e} (<= {eps_dual:.3e})  rho {:.3e}",
                ws.rho
            );
        }

        if r_prim <= eps_prim && r_dual <= eps_dual {
            status = Status::Solved;
            iterations = iteration;
            break;
        }
        if last {
            break;
        }

        // Penalty update. Both ratios are relative, so the estimate is scale free; a
        // vanishing denominator means that residual is already at rounding level and
        // carries no information about which way to move.
        let rel_prim = if scale_prim > residual_floor {
            r_prim / scale_prim
        } else {
            0.0
        };
        let rel_dual = if scale_dual > residual_floor {
            r_dual / scale_dual
        } else {
            0.0
        };
        if rel_dual > 0.0 && rel_prim > 0.0 {
            let candidate = (ws.rho * (rel_prim / rel_dual).sqrt()).clamp(RHO_MIN, rho_limit);
            let drift = candidate / ws.rho;
            if drift > ADAPTIVE_RHO_TOLERANCE || drift < 1.0 / ADAPTIVE_RHO_TOLERANCE {
                ws.rho = candidate;
                refresh_factor(
                    a,
                    lambda_max_bound,
                    settings.eps_abs,
                    ws,
                    factor,
                    &mut refactorizations,
                );
            }
        }
    }

    Outcome {
        status,
        iterations,
        refactorizations,
        primal_residual,
        dual_residual,
        rho: ws.rho,
        used_f32: ws.factor_used_f32,
    }
}

/// Rebuild the numeric factorization if `rho` moved or `A`'s values changed.
///
/// The symbolic phase is untouched -- only `K`'s values depend on `rho` and `coef`, and
/// its pattern is fixed for the whole fit.
fn refresh_factor(
    a: &Constraints<'_>,
    lambda_max_bound: f64,
    eps_abs: f64,
    ws: &mut Workspace,
    factor: &mut KktFactor,
    refactorizations: &mut u32,
) {
    if ws.factor_rho == ws.rho {
        return;
    }
    // `P = I` floors the spectrum at `1 + sigma`, so this is the exact condition number,
    // not an estimate, and the f32 forward error it implies is a measured linear function
    // of it. Storing the factor in f32 halves the traffic of the triangular solves, which
    // dominate the per-iteration cost once the factorization is amortized over a threshold
    // sequence -- but only while the resulting error stays well inside `eps_abs`.
    let condition = 1.0 + ws.rho * lambda_max_bound / (1.0 + SIGMA);
    let use_f32 = F32_ERROR_PER_COND * condition <= F32_ERROR_MARGIN * eps_abs;
    factor.refactor(SIGMA, ws.rho, a.edges, a.coef, use_f32);
    ws.factor_rho = ws.rho;
    ws.factor_used_f32 = use_f32;
    *refactorizations += 1;
}

/// The largest `rho` whose linear system f64 can still solve accurately enough to reach
/// `eps_abs`.
///
/// `rho` trades primal feasibility against dual optimality, and the usual heuristic will
/// happily drive it to its bound. But `rho` also sets the conditioning of the system being
/// solved -- `cond(K) = 1 + rho lambda_max(A^T A) / (1 + sigma)`, exactly, because `P = I`
/// floors the spectrum -- and a solve is only accurate to about `eps * cond(K)`. Past the
/// point where that product reaches the requested tolerance, raising `rho` stops buying
/// primal feasibility and merely installs a floor under the dual residual, so the solve
/// runs to its iteration limit having already found the answer. Capping `rho` is what
/// keeps the termination test reachable.
fn rho_ceiling(lambda_max_bound: f64, eps_abs: f64) -> f64 {
    if lambda_max_bound <= 0.0 || !lambda_max_bound.is_finite() {
        return RHO_MAX;
    }
    let condition_limit = F32_ERROR_MARGIN * eps_abs / F64_ERROR_PER_COND;
    let ceiling = (condition_limit - 1.0) * (1.0 + SIGMA) / lambda_max_bound;
    ceiling.clamp(RHO_MIN, RHO_MAX)
}
