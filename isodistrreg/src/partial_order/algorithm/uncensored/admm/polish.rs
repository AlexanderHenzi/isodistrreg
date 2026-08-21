//! Closed-form polishing of an ADMM iterate.
//!
//! An active constraint row `(i, j)` with coefficients `(a, b)` says `a u_i + b u_j = 0`,
//! which pins the *ratio* of the two variables. The active rows therefore partition the
//! covariates into components on which the solution is determined up to a single scalar:
//! fix `g` at one node of a component and propagate it along the active edges, and the
//! component's contribution to the objective collapses to a scalar quadratic in `t`,
//! minimized at
//!
//! ```text
//!     t* = - (sum_{c in C} q_c g_c) / (sum_{c in C} g_c^2),      u_c = t* g_c
//! ```
//!
//! For stochastic dominance every row is `+-s (x_j - x_i)`, so `g` is proportional to
//! `sqrt_weight` and this is the weight-weighted block mean -- the value PAVA would pool
//! to. A singleton component gives `t* = -q_c`, the unconstrained minimizer, as it must.
//!
//! Since `P = I` is positive definite, the true optimum *is* the minimizer over the
//! equalities induced by the correct active set, so identifying that set exactly makes
//! this step exact. ADMM identifies it structurally: `z = max(0, w)` and the entries it
//! clips are the active rows, recorded bit for bit in `Workspace::clipped`.
//!
//! Two safeguards make the step never-worse rather than a gamble:
//!
//!  * a **repair loop** -- any row the polished point violates is added to the active set
//!    and the polish repeated. The active set only grows, so this terminates, and it is
//!    exactly PAVA's pooling of an out-of-order block.
//!  * an **accept gate** -- the result replaces the iterate only if it is feasible and has
//!    a strictly lower objective. For two feasible points,
//!    `||u - u*|| <= sqrt(2 (f(u) - f*))`, so a point that passes is provably no further
//!    from the optimum than the one it replaces.
//!
//! [`Polisher::peel_multipliers`] closes the loop on the dual side: the spanning forest
//! the ratio propagation walks makes the stationarity system triangular, so the polished
//! point's exact multipliers come out in `O(n + m)`. Nonnegative multipliers plus a
//! rounding-level duality gap certify the polished point as *the* optimum; a negative one
//! names the row the active set should drop. The solver builds its exact-finish
//! predictor/corrector out of exactly these two signals.
use super::constraints::Constraints;

/// Feasibility slack allowed of the polished point, absolute in energy coordinates.
const FEASIBILITY_TOLERANCE: f64 = 1e-9;
/// Relative slack when checking that a cycle's propagated ratios agree.
const CONSISTENCY_TOLERANCE: f64 = 1e-9;

/// Scratch for polishing, sized once per fit and reused across thresholds.
pub(crate) struct Polisher {
    /// The candidate solution; valid after `polish` returns `true`.
    pub(crate) u: Vec<f64>,
    /// Multipliers on the active spanning forest; valid after `peel_multipliers`.
    pub(crate) multipliers: Vec<f64>,
    active: Vec<bool>,
    /// Per-component scalar `g`, and the component id each covariate belongs to.
    g: Vec<f64>,
    component: Vec<usize>,
    consistent: Vec<bool>,
    numerator: Vec<f64>,
    denominator: Vec<f64>,
    /// Whether every `q` entry of the component is exactly zero.
    q_zero: Vec<bool>,
    stack: Vec<usize>,
    adj_ptr: Vec<usize>,
    adj_cursor: Vec<usize>,
    adj: Vec<usize>,
    au: Vec<f64>,
    /// The active edge over which each covariate was first reached in `propagate_ratios`,
    /// `usize::MAX` at component roots -- a spanning forest of the active subgraph.
    parent_edge: Vec<usize>,
    /// Covariates in discovery order; every node appears after its forest parent, so the
    /// reverse order peels leaves first.
    discovery: Vec<usize>,
    /// Stationarity residual scratch for the peel.
    residual: Vec<f64>,
}

/// What `Polisher::peel_multipliers` found.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Peel {
    /// How many multipliers were strictly negative (clamped to zero in the stored vector).
    pub(crate) negative: usize,
    /// The row with the most negative multiplier, `usize::MAX` when none was negative.
    pub(crate) worst_row: usize,
}

/// What the polish did, for diagnostics and for the caller's acceptance logic.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct PolishReport {
    /// Feasible and strictly better than `raw_u` -- the never-worse replacement gate.
    pub(crate) accepted: bool,
    /// Feasible (within tolerance) regardless of objective. A feasible candidate can
    /// still be certified optimal by a duality gap even when it does not beat `raw_u`,
    /// e.g. when `raw_u` is slightly infeasible with a lower objective.
    pub(crate) feasible: bool,
    pub(crate) repair_rounds: u32,
}

impl Polisher {
    /// The active set the last `polish` call ended on, after any repair rounds.
    pub(crate) fn active(&self) -> &[bool] {
        &self.active
    }

    /// Whether the last polish left any component with contradictory ratios.
    pub(crate) fn any_inconsistent(&self) -> bool {
        self.consistent.iter().any(|&ok| !ok)
    }

    pub(crate) fn new(n: usize, m: usize) -> Self {
        Self {
            u: vec![0.0; n],
            multipliers: vec![0.0; m],
            active: vec![false; m],
            g: vec![0.0; n],
            component: vec![usize::MAX; n],
            consistent: Vec::with_capacity(n),
            numerator: Vec::with_capacity(n),
            denominator: Vec::with_capacity(n),
            q_zero: Vec::with_capacity(n),
            stack: Vec::with_capacity(n),
            adj_ptr: vec![0; n + 1],
            adj_cursor: vec![0; n + 1],
            adj: vec![0; 2 * m],
            au: vec![0.0; m],
            parent_edge: vec![usize::MAX; n],
            discovery: Vec::with_capacity(n),
            residual: vec![0.0; n],
        }
    }

    /// Build the adjacency of the currently active rows, as CSR over covariates.
    fn build_adjacency(&mut self, a: &Constraints<'_>) {
        let n = a.n;
        self.adj_ptr[..=n].fill(0);
        for (r, &(i, j)) in a.edges.iter().enumerate() {
            if self.active[r] {
                self.adj_ptr[i] += 1;
                self.adj_ptr[j] += 1;
            }
        }
        let mut running = 0usize;
        for slot in self.adj_ptr[..n].iter_mut() {
            let count = *slot;
            *slot = running;
            running += count;
        }
        self.adj_ptr[n] = running;
        self.adj_cursor[..=n].copy_from_slice(&self.adj_ptr[..=n]);
        for (r, &(i, j)) in a.edges.iter().enumerate() {
            if self.active[r] {
                self.adj[self.adj_cursor[i]] = r;
                self.adj_cursor[i] += 1;
                self.adj[self.adj_cursor[j]] = r;
                self.adj_cursor[j] += 1;
            }
        }
    }

    /// Assign every covariate a component and a ratio `g`, flagging components whose
    /// active edges pin contradictory ratios around a cycle.
    fn propagate_ratios(&mut self, a: &Constraints<'_>) {
        let n = a.n;
        self.component[..n].fill(usize::MAX);
        self.consistent.clear();
        self.parent_edge[..n].fill(usize::MAX);
        self.discovery.clear();
        let mut n_components = 0usize;

        for root in 0..n {
            if self.component[root] != usize::MAX {
                continue;
            }
            let id = n_components;
            n_components += 1;
            self.consistent.push(true);
            self.component[root] = id;
            self.g[root] = 1.0;
            self.discovery.push(root);
            self.stack.clear();
            self.stack.push(root);

            while let Some(c) = self.stack.pop() {
                for slot in self.adj_ptr[c]..self.adj_ptr[c + 1] {
                    let r = self.adj[slot];
                    let (i, j) = a.edges[r];
                    let (coef_i, coef_j) = a.coef[r];
                    let (other, coef_c, coef_other) = if c == i {
                        (j, coef_i, coef_j)
                    } else {
                        (i, coef_j, coef_i)
                    };

                    // a_c g_c + a_other g_other = 0. Both coefficients are nonzero for
                    // positive weights; a zero or non-finite one means the ratio is not
                    // determined, so the component is not polishable.
                    if coef_other == 0.0 || !coef_other.is_finite() || !coef_c.is_finite() {
                        self.consistent[id] = false;
                        continue;
                    }
                    let implied = -coef_c * self.g[c] / coef_other;
                    if !implied.is_finite() {
                        self.consistent[id] = false;
                        continue;
                    }
                    if self.component[other] == usize::MAX {
                        self.component[other] = id;
                        self.g[other] = implied;
                        self.parent_edge[other] = r;
                        self.discovery.push(other);
                        self.stack.push(other);
                    } else {
                        let residual = coef_c * self.g[c] + coef_other * self.g[other];
                        let scale = (coef_c * self.g[c])
                            .abs()
                            .max((coef_other * self.g[other]).abs());
                        if residual.abs() > CONSISTENCY_TOLERANCE * scale.max(1.0) {
                            self.consistent[id] = false;
                        }
                    }
                }
            }
        }

        self.numerator.clear();
        self.numerator.resize(n_components, 0.0);
        self.denominator.clear();
        self.denominator.resize(n_components, 0.0);
    }

    /// Minimize the objective over the ray each component is confined to.
    fn solve_components(&mut self, q: &[f64], raw_u: &[f64]) {
        for slot in self.numerator.iter_mut() {
            *slot = 0.0;
        }
        for slot in self.denominator.iter_mut() {
            *slot = 0.0;
        }
        self.q_zero.clear();
        self.q_zero.resize(self.numerator.len(), true);
        for (c, (&q_c, &g_c)) in q.iter().zip(&self.g).enumerate() {
            let id = self.component[c];
            self.numerator[id] += q_c * g_c;
            self.denominator[id] += g_c * g_c;
            self.q_zero[id] &= q_c == 0.0;
        }
        for (c, u_c) in self.u.iter_mut().enumerate() {
            let id = self.component[c];
            let denominator = self.denominator[id];
            if !self.consistent[id] || denominator <= 0.0 || !denominator.is_finite() {
                // Contradictory ratios around a cycle admit only `u = 0` on the
                // component, so zero *is* the exact minimizer over these equalities. When
                // the component's data is all zero as well -- the hazard-rate path's
                // fully-passed region, where the fitted survival is exactly zero and the
                // ratio rows have degenerated -- that point is also stationary with zero
                // multipliers, so a certificate over it holds. With nonzero data the zero
                // point would not be stationary and could sit far from the iterate, so
                // the solver's own value is the safer fallback there; the acceptance
                // gates judge either choice.
                // (peel_multipliers skips these components, so their multipliers stay
                // zero and any certificate computed from them stays honest.)
                *u_c = if self.q_zero[id] { 0.0 } else { raw_u[c] };
                continue;
            }
            let t = -self.numerator[id] / denominator;
            *u_c = if t.is_finite() {
                t * self.g[c]
            } else {
                raw_u[c]
            };
        }
    }

    /// Exact multipliers of the polished point, by peeling leaves of the active forest.
    ///
    /// Stationarity of the QP is `u + q = A^T lambda` with `lambda` supported on the
    /// active rows. Restricted to the spanning forest recorded by `propagate_ratios` that
    /// system is triangular: reverse discovery order visits every node after all of its
    /// forest children, so each node's remaining residual determines its parent edge's
    /// multiplier in O(1). At a component root the residual cancels to rounding, because
    /// the pooled value `t*` is exactly the zero of `sum_c g_c (u_c + q_c)`. Redundant
    /// active edges -- the ones that close cycles -- keep multiplier zero, which is valid:
    /// dropping them from the support changes neither the partition nor the pooled values.
    ///
    /// Returns how many multipliers came out strictly negative and which row's was the
    /// most negative (`usize::MAX` when none); negatives are clamped to zero in
    /// `multipliers`, so the stored vector is dual-feasible either way and a duality gap
    /// computed from it is a genuine bound. A negative multiplier is the dual's way of
    /// saying that edge should not be active (a pool that wants to split); dropping the
    /// most negative one and re-solving is the classical primal active-set step, which --
    /// unlike dropping all of them at once -- cannot cycle.
    ///
    /// Components flagged inconsistent are skipped; their multipliers stay zero, their
    /// stationarity residual stays in the gap, and a certificate over them simply fails.
    ///
    /// `blocked` rows are excluded from the `worst_row` selection (they still count as
    /// negative): the caller uses this to never drop the same row twice within one finish
    /// attempt, which is what keeps the exchange from cycling through a degenerate set of
    /// splits.
    pub(crate) fn peel_multipliers(
        &mut self,
        a: &Constraints<'_>,
        q: &[f64],
        blocked: Option<&[bool]>,
    ) -> Peel {
        let n = a.n;
        self.multipliers[..a.m()].fill(0.0);
        for (residual_c, (&u_c, &q_c)) in self.residual[..n].iter_mut().zip(self.u.iter().zip(q)) {
            *residual_c = u_c + q_c;
        }

        let mut peel = Peel {
            negative: 0,
            worst_row: usize::MAX,
        };
        let mut worst = 0.0f64;
        for &c in self.discovery.iter().rev() {
            let r = self.parent_edge[c];
            if r == usize::MAX || !self.consistent[self.component[c]] {
                continue;
            }
            let (i, j) = a.edges[r];
            let (coef_i, coef_j) = a.coef[r];
            // `c` was discovered over `r`, so its coefficient passed the nonzero check in
            // `propagate_ratios`; the division is safe.
            let (coef_c, other, coef_other) = if c == i {
                (coef_i, j, coef_j)
            } else {
                (coef_j, i, coef_i)
            };
            let lambda = self.residual[c] / coef_c;
            // Propagate the *unclamped* value so the remaining tree multipliers stay the
            // exact solution of the forest system; only the stored vector is clamped.
            self.residual[other] -= coef_other * lambda;
            if lambda > 0.0 {
                self.multipliers[r] = lambda;
            } else if lambda < 0.0 {
                peel.negative += 1;
                if lambda < worst && blocked.is_none_or(|blocked| !blocked[r]) {
                    worst = lambda;
                    peel.worst_row = r;
                }
            }
        }
        peel
    }
}

fn objective(u: &[f64], q: &[f64]) -> f64 {
    let mut quadratic = 0.0;
    let mut linear = 0.0;
    for (&u_c, &q_c) in u.iter().zip(q) {
        quadratic += u_c * u_c;
        linear += q_c * u_c;
    }
    0.5 * quadratic + linear
}

/// Try to replace `raw_u` with the exact minimizer over the active set.
///
/// `active_seed(r)` supplies the starting active set -- the projection clips recorded in
/// `Workspace::clipped`, or the support of a dual iterate's multipliers. Returns a
/// report; on `feasible` the candidate point is in `p.u`. Otherwise `p.u` is left in an
/// unspecified state and the caller keeps `raw_u`.
pub(crate) fn polish(
    a: &Constraints<'_>,
    q: &[f64],
    active_seed: impl Fn(usize) -> bool,
    raw_u: &[f64],
    p: &mut Polisher,
) -> PolishReport {
    let m = a.m();
    for (r, active_r) in p.active[..m].iter_mut().enumerate() {
        *active_r = active_seed(r);
    }

    let mut report = PolishReport::default();
    // The active set only grows, so each round either finishes or adds a row; `m + 1`
    // rounds is therefore a bound, not a heuristic cutoff.
    for round in 0..=m {
        p.build_adjacency(a);
        p.propagate_ratios(a);
        p.solve_components(q, raw_u);

        a.mul(&p.u, &mut p.au);
        let mut violated = false;
        for (r, &au_r) in p.au.iter().enumerate() {
            if au_r < -FEASIBILITY_TOLERANCE && !p.active[r] {
                p.active[r] = true;
                violated = true;
            }
        }
        if !violated {
            report.repair_rounds = round as u32;
            break;
        }
    }

    // Accept only a feasible point with a strictly better objective; otherwise the
    // solver's own iterate stands.
    report.feasible = p.au.iter().all(|&v| v >= -FEASIBILITY_TOLERANCE);
    report.accepted = report.feasible && objective(&p.u, q) < objective(raw_u, q);
    report
}
