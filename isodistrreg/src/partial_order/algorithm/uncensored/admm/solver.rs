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
//!
//! The iteration is wrapped in certified exact-finish attempts (see [`try_finish`]): a
//! predictor before the first iteration and correctors at failing residual checks, each
//! accepting only on a duality gap of at most `0.5 * eps_abs^2` -- by strong convexity a
//! guarantee of `||u - u*|| <= eps_abs`, strictly stronger than the residual test. Warm
//! starts make consecutive thresholds share their active structure, so most solves never
//! reach the first iteration at all, and the ones that do stop as soon as the projection
//! clips settle rather than when the residuals finish crawling.
use super::constraints::Constraints;
use super::kkt::KktFactor;
use super::polish::{self, Polisher};

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

/// Active-set refinement rounds for the predictor attempt at a solve's entry.
///
/// Each round polishes a candidate active set, peels its exact multipliers, and either
/// certifies or drops the most negative multiplier's row -- the dual's exact signal for a
/// pool that wants to split, applied one row at a time exactly as the classical primal
/// active-set method does (dropping several at once was observed to limit-cycle).
/// Consecutive thresholds differ by a bounded cascade of such splits plus the merges the
/// repair loop restores, so this cap either finishes the threshold outright or
/// establishes that the structure genuinely moved and the splitting iteration should run.
const PREDICTOR_ROUNDS: u32 = 16;

/// Refinement rounds for corrector attempts at failed residual checks.
///
/// Mid-iteration the clip record is noisy, and measured across instances a corrector that
/// cannot certify within a few exchanges almost never certifies with more -- the deep
/// marches are the predictor's job -- so the in-loop probes stay shallow.
const CORRECTOR_ROUNDS: u32 = 4;

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
    /// Primal residual `||A u - z||_inf` at the returned iterate. Zero when the solve
    /// finished through a certified exact polish, whose acceptance implies feasibility.
    pub(crate) primal_residual: f64,
    /// Dual residual `||u + q + A^T y||_inf` at the returned iterate. Zero when the solve
    /// finished through a certified exact polish.
    pub(crate) dual_residual: f64,
    /// The penalty the solve ended on, carried into the next threshold.
    pub(crate) rho: f64,
    /// Whether the final factorization was stored in f32.
    pub(crate) used_f32: bool,
    /// The iterate in `Workspace::x` is a polished point whose duality gap certifies
    /// `||u - u*|| <= eps_abs`; the caller's own polish pass would be redundant.
    pub(crate) finished_exact: bool,
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
    /// Refined active-set candidate handed between finish rounds; lives here rather than
    /// in the `Polisher` because the polish call mutably borrows the polisher while the
    /// seed is read.
    seed_scratch: Vec<bool>,
    /// Rows already dropped in the current finish attempt; never dropped twice, which
    /// keeps the active-set exchange from cycling through a degenerate set of splits.
    dropped_scratch: Vec<bool>,
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
            seed_scratch: vec![false; m],
            dropped_scratch: vec![false; m],
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

/// `f(u) - g(lambda)`: the duality gap of a primal point against explicit multipliers,
/// with `g(lambda) = -0.5 ||A^T lambda - q||^2` the dual of `min 0.5||u||^2 + q^T u`
/// s.t. `A u >= 0`. Valid as an optimality bound for any `lambda >= 0` by weak duality.
///
/// Summed with Kahan compensation: acceptance compares the gap against
/// `0.5 * eps_abs^2 ~ 5e-11` while `f` and `g` are sums of `n` order-one terms, so a
/// plain sequential sum would drown the threshold in its own rounding for a few thousand
/// covariates. The per-coordinate terms `0.5 u_c^2 + q_c u_c + 0.5 v_c^2` are formed
/// first -- each is order one, and the near-total cancellation happens across
/// coordinates, exactly where the compensation acts.
///
/// `v_scratch` is an `n`-length buffer; its contents on entry are ignored.
fn certified_gap(
    a: &Constraints<'_>,
    q: &[f64],
    u: &[f64],
    multipliers: &[f64],
    v_scratch: &mut [f64],
) -> f64 {
    a.mul_transpose(multipliers, v_scratch);
    let mut sum = 0.0f64;
    let mut compensation = 0.0f64;
    for ((&u_c, &q_c), &aty_c) in u.iter().zip(q).zip(v_scratch.iter()) {
        let v_c = aty_c - q_c;
        let term = 0.5 * u_c * u_c + q_c * u_c + 0.5 * v_c * v_c;
        let adjusted = term - compensation;
        let next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    sum
}

/// Try to finish the solve exactly: polish a candidate active set and accept it only on a
/// rigorous duality-gap certificate.
///
/// The first candidate is the projection-clip record (the previous threshold's active set
/// on a fresh solve). Polishing it yields the exact minimizer over its equalities, and
/// peeling the active forest yields that point's exact multipliers, so when the set is
/// right the gap is rounding-level and certifies immediately -- no iteration, and on the
/// predictor call no factorization either. When multipliers come out negative, their
/// support is the dual's exact instruction for which rows to drop (a pool splitting as
/// the threshold moves), and the next round polishes that; the repair loop inside
/// `polish` re-adds any row the split point violates. Either the process closes the gap
/// within a few rounds or the structure genuinely moved and the splitting iteration runs.
///
/// Acceptance at `gap <= 0.5 * eps_abs^2` converts, by strong convexity of the
/// unit-quadratic objective, into `||u - u*|| <= eps_abs` -- a guarantee at least as
/// strong as the residual test, typically available long before the residuals crawl
/// under their tolerances.
///
/// On success `ws.x` holds the certified point and `ws.clipped` its active set (the seed
/// for the next threshold's attempt); on failure both are untouched.
fn try_finish(
    a: &Constraints<'_>,
    q: &[f64],
    ws: &mut Workspace,
    polisher: &mut Polisher,
    gap_tolerance: f64,
    rounds: u32,
) -> Option<f64> {
    let m = a.m();
    ws.dropped_scratch[..m].fill(false);
    for round in 0..rounds {
        let report = if round == 0 {
            polish::polish(a, q, |r| ws.clipped[r], &ws.x, polisher)
        } else {
            polish::polish(a, q, |r| ws.seed_scratch[r], &ws.x, polisher)
        };
        let peel = polisher.peel_multipliers(a, q, Some(&ws.dropped_scratch));
        if report.feasible {
            let gap = certified_gap(a, q, &polisher.u, &polisher.multipliers, &mut ws.aty);
            if gap <= gap_tolerance {
                ws.x.copy_from_slice(&polisher.u);
                ws.clipped.copy_from_slice(polisher.active());
                return Some(gap);
            }
            if peel.negative == 0 {
                // Stationary with nonnegative multipliers yet a large gap: some component
                // was not polishable (inconsistent ratios with data on it), which more
                // rounds cannot fix.
                return None;
            }
        } else if peel.negative == 0 && !polisher.any_inconsistent() {
            // Infeasible with no refinement signal: nothing left to drop.
            return None;
        }
        if polisher.any_inconsistent() {
            // The clip record can glue the degenerate zero-survival region to live
            // covariates. Keep exactly the edges the peel priced strictly positive: that
            // disintegrates every inconsistent component wholesale and lets the repair
            // loop rebuild the region from consistent equalities.
            for (seed_r, &lambda_r) in ws.seed_scratch[..m].iter_mut().zip(&polisher.multipliers) {
                *seed_r = lambda_r > 0.0;
            }
        } else {
            // The classical active-set step: drop only the most negative multiplier's
            // row. Dropping several at once re-merges via the repair loop and was
            // observed to limit-cycle; one at a time makes each round a strict exchange,
            // and never re-dropping a row the repair loop re-added keeps the exchange
            // out of degenerate split cycles.
            if peel.worst_row == usize::MAX {
                // Every negative multiplier's row was already dropped once: no fresh
                // exchange is left to try.
                return None;
            }
            for (seed_r, &active_r) in ws.seed_scratch[..m].iter_mut().zip(polisher.active()) {
                *seed_r = active_r;
            }
            ws.seed_scratch[peel.worst_row] = false;
            ws.dropped_scratch[peel.worst_row] = true;
        }
    }
    None
}

/// Solve `min 0.5 u^T u + q^T u  s.t.  A u >= 0`, leaving the primal iterate in `ws.x`.
///
/// `ws` carries the warm start in and the solution out. `factor` must have been built
/// with `KktFactor::new` for the same `n` and edge list; this routine refreshes its
/// numeric values whenever `rho` moves or `Workspace::invalidate_factor` was called.
///
/// `polisher` drives the exact-finish attempts: one before the first iteration
/// (consecutive thresholds usually share their active set, so the previous solution's
/// structure plus the new `q` is often already the answer -- in which case the
/// factorization is never touched at all), and one at every residual check that fails.
pub(crate) fn solve(
    a: &Constraints<'_>,
    q: &[f64],
    settings: &Settings,
    ws: &mut Workspace,
    factor: &mut KktFactor,
    polisher: &mut Polisher,
) -> Outcome {
    debug_assert_eq!(q.len(), a.n);
    let gap_tolerance = 0.5 * settings.eps_abs * settings.eps_abs;

    // Predictor: the warm-started structure certifies without a single iteration on most
    // thresholds. Runs before the spectral bound and the factorization refresh, so an
    // exact finish here skips the numeric factorization entirely -- on the hazard-rate
    // path, where `A` changes every threshold, that is the dominant saving.
    if let Some(gap) = try_finish(a, q, ws, polisher, gap_tolerance, PREDICTOR_ROUNDS) {
        if settings.verbose {
            eprintln!("  admm predictor finish, certified gap {gap:.3e}");
        }
        return Outcome {
            status: Status::Solved,
            iterations: 0,
            refactorizations: 0,
            primal_residual: 0.0,
            dual_residual: 0.0,
            rho: ws.rho,
            used_f32: ws.factor_used_f32,
            finished_exact: true,
        };
    }

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

    // Exponential backoff for corrector attempts. When the iterate's structure is still
    // far from the answer -- degenerate near-equal pools mid-grind -- every attempt fails
    // the certificate, and paying several edge-list passes per residual check adds up.
    // Doubling the wait after each failure keeps a hard solve's attempt cost logarithmic
    // in its iteration count while an easy solve still exits on the first check.
    let mut checks_failed = 0u32;
    let mut next_attempt = 1u32;
    let mut previous_worst_ratio = f64::INFINITY;

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

        // Corrector: the residuals are still crawling, but the *structure* of the iterate
        // -- which rows the projection clips -- typically settles orders of magnitude
        // earlier, and once it has, the polished point is exact and certifiable. This is
        // what cuts off the slow tail of the splitting iteration.
        checks_failed += 1;
        if checks_failed >= next_attempt || last {
            next_attempt = checks_failed
                .saturating_mul(2)
                .max(next_attempt.saturating_mul(2));
            if let Some(gap) = try_finish(a, q, ws, polisher, gap_tolerance, CORRECTOR_ROUNDS) {
                if settings.verbose {
                    eprintln!(
                        "  admm iter {iteration:>5}  corrector finish, certified gap {gap:.3e}"
                    );
                }
                return Outcome {
                    status: Status::Solved,
                    iterations: iteration,
                    refactorizations,
                    primal_residual: 0.0,
                    dual_residual: 0.0,
                    rho: ws.rho,
                    used_f32: ws.factor_used_f32,
                    finished_exact: true,
                };
            }
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
        // An endgame plateau overrides the adoption deadband. The deadband exists to
        // avoid refactorization churn while the iteration is making progress -- but a
        // residual ratio stuck at, say, 5 sits exactly inside it, and the iteration then
        // decays at a fraction of a percent per sweep until the budget burns. Barely
        // improving residuals mean the current `rho` has nothing more to give. The
        // near-convergence gate matters: early in a solve the residuals wander without
        // shrinking monotonically, and treating that as a plateau chases `rho` far past
        // any useful value.
        let worst_ratio =
            (r_prim / eps_prim.max(residual_floor)).max(r_dual / eps_dual.max(residual_floor));
        let plateaued = worst_ratio > previous_worst_ratio * 0.95 && worst_ratio <= 32.0;
        previous_worst_ratio = worst_ratio;
        if rel_dual > 0.0 && rel_prim > 0.0 {
            let candidate = (ws.rho * (rel_prim / rel_dual).sqrt()).clamp(RHO_MIN, rho_limit);
            let drift = candidate / ws.rho;
            let outside_deadband =
                drift > ADAPTIVE_RHO_TOLERANCE || drift < 1.0 / ADAPTIVE_RHO_TOLERANCE;
            let meaningful = drift > 1.05 || drift < 1.0 / 1.05;
            if outside_deadband || (plateaued && meaningful) {
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
        finished_exact: false,
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
