//! Hildreth dual coordinate ascent: a rigorous accuracy certificate, and a fallback.
//!
//! The Lagrangian of `min 0.5 u^T u + q^T u  s.t.  A u >= 0` is stationary at
//! `u = A^T lambda - q`, and substituting it back leaves the dual
//!
//! ```text
//!     max_{lambda >= 0}  g(lambda) = -0.5 || A^T lambda - q ||^2
//! ```
//!
//! Weak duality makes `g(lambda)` a lower bound on the primal optimum for *any*
//! `lambda >= 0`, with no convergence assumption. Paired with a feasible primal point `u`
//! that gives a certified gap, and strong convexity turns the gap into a distance:
//! `||u - u*|| <= sqrt(2 (f(u) - g(lambda)))`. That is the only way this solver can say
//! how wrong a non-converged answer is, rather than guess.
//!
//! Hildreth's method maximizes `g` one multiplier at a time. Along coordinate `r` the
//! exact step is `-a_r^T v / ||a_r||^2` followed by a clamp at zero, where `v` is the
//! running residual `A^T lambda - q`. Because every row of `A` has exactly two nonzeros,
//! both the step and the residual update are O(1), so a full sweep is O(m) with no
//! factorization at all -- which is why it still makes progress where the ADMM iteration
//! has stalled.
use super::constraints::Constraints;

/// Feasibility slack allowed of a candidate point, absolute in energy coordinates.
const FEASIBILITY_TOLERANCE: f64 = 1e-9;

fn objective(u: &[f64], q: &[f64]) -> f64 {
    0.5 * u.iter().map(|value| value * value).sum::<f64>()
        + u.iter()
            .zip(q)
            .map(|(value, cost)| value * cost)
            .sum::<f64>()
}

/// Coordinate sweeps to spend when the splitting iteration has stalled.
///
/// A sweep is `O(m)` with no factorization, so this is cheap next to the iterations that
/// already failed; it only runs on a solve that did not converge.
pub(crate) const STALL_SWEEPS: u32 = 200;

/// Dual iterate and the running primal residual it implies.
pub(crate) struct Hildreth {
    lambda: Vec<f64>,
    /// `A^T lambda - q`, which is exactly the primal point the multipliers imply.
    v: Vec<f64>,
    row_norm_squared: Vec<f64>,
    row_scratch: Vec<f64>,
}

impl Hildreth {
    pub(crate) fn new(n: usize, m: usize) -> Self {
        Self {
            lambda: vec![0.0; m],
            v: vec![0.0; n],
            row_norm_squared: vec![0.0; m],
            row_scratch: vec![0.0; m],
        }
    }

    /// Start from `lambda = max(0, -y)`, the multipliers implied by an ADMM dual iterate.
    ///
    /// The ADMM stationarity condition is `u + q + A^T y = 0` while this dual uses
    /// `u + q - A^T lambda = 0`, so the two differ by a sign.
    pub(crate) fn seed_from_admm_dual(&mut self, a: &Constraints<'_>, q: &[f64], y: &[f64]) {
        for (lambda_r, &y_r) in self.lambda.iter_mut().zip(y) {
            *lambda_r = if y_r < 0.0 { -y_r } else { 0.0 };
        }
        self.refresh(a, q);
    }

    /// Start from `lambda = 0`.
    #[cfg(test)]
    pub(crate) fn seed_zero(&mut self, a: &Constraints<'_>, q: &[f64]) {
        self.lambda.fill(0.0);
        self.refresh(a, q);
    }

    fn refresh(&mut self, a: &Constraints<'_>, q: &[f64]) {
        a.mul_transpose(&self.lambda, &mut self.v);
        for (v_c, &q_c) in self.v.iter_mut().zip(q) {
            *v_c -= q_c;
        }
        for (norm_r, &(coef_i, coef_j)) in self.row_norm_squared.iter_mut().zip(a.coef) {
            *norm_r = coef_i * coef_i + coef_j * coef_j;
        }
    }

    /// Run `sweeps` passes of coordinate ascent. Monotone in the dual objective.
    pub(crate) fn run(&mut self, a: &Constraints<'_>, sweeps: u32) {
        for _ in 0..sweeps {
            for (r, (&(i, j), &(coef_i, coef_j))) in a.edges.iter().zip(a.coef).enumerate() {
                let norm = self.row_norm_squared[r];
                if norm <= 0.0 || !norm.is_finite() {
                    continue;
                }
                let gradient = coef_i * self.v[i] + coef_j * self.v[j];
                let candidate = self.lambda[r] - gradient / norm;
                let next = if candidate > 0.0 { candidate } else { 0.0 };
                let step = next - self.lambda[r];
                if step == 0.0 {
                    continue;
                }
                self.lambda[r] = next;
                self.v[i] += coef_i * step;
                self.v[j] += coef_j * step;
            }
        }
    }

    /// `g(lambda) = -0.5 ||A^T lambda - q||^2`: a lower bound on the primal optimum,
    /// valid at any iterate.
    pub(crate) fn dual_bound(&self) -> f64 {
        -0.5 * self.v.iter().map(|value| value * value).sum::<f64>()
    }

    /// The primal point the multipliers imply, `u = A^T lambda - q`.
    #[cfg(test)]
    pub(crate) fn primal(&self) -> &[f64] {
        &self.v
    }

    #[cfg(test)]
    pub(crate) fn multipliers(&self) -> &[f64] {
        &self.lambda
    }

    /// A rigorous upper bound on this point's suboptimality, `f(u) - g(lambda)`.
    ///
    /// Weak duality makes it valid at any dual-feasible `lambda`, so it reports how far a
    /// non-converged answer actually is rather than assuming. Meaningful only for a
    /// primal-feasible `u`; an infeasible one can undercut the bound and return a negative
    /// number.
    pub(crate) fn certified_gap(&self, q: &[f64], u: &[f64]) -> f64 {
        objective(u, q) - self.dual_bound()
    }

    /// Overwrite `u` with the point these multipliers imply, but only if that point is
    /// feasible and strictly better.
    ///
    /// Two feasible points obey `||u - u*|| <= sqrt(2 (f(u) - f*))`, so accepting a
    /// lower-objective feasible point can only move `u` closer to the optimum. Returns
    /// whether the swap happened.
    pub(crate) fn take_if_better(&mut self, a: &Constraints<'_>, q: &[f64], u: &mut [f64]) -> bool {
        a.mul(&self.v, &mut self.row_scratch);
        let feasible = self
            .row_scratch
            .iter()
            .all(|&row| row >= -FEASIBILITY_TOLERANCE);
        // Written as a positive comparison bound to a name: a NaN objective must fail the
        // test rather than pass it, which `>=` would not do.
        let improves = objective(&self.v, q) < objective(u, q);
        if !feasible || !improves {
            return false;
        }
        u.copy_from_slice(&self.v);
        true
    }
}

#[cfg(test)]
mod test {
    use super::Hildreth;
    use crate::partial_order::algorithm::uncensored::admm::constraints::Constraints;

    fn objective(u: &[f64], q: &[f64]) -> f64 {
        0.5 * u.iter().map(|v| v * v).sum::<f64>()
            + u.iter().zip(q).map(|(a, b)| a * b).sum::<f64>()
    }

    /// With no constraints the dual is maximized at `lambda = 0` and the bound is exactly
    /// the unconstrained optimum `-0.5 ||q||^2`.
    #[test]
    fn unconstrained_bound_is_tight() {
        let a = Constraints {
            edges: &[],
            coef: &[],
            n: 3,
        };
        let q = [1.0, -2.0, 0.5];
        let mut h = Hildreth::new(3, 0);
        h.seed_zero(&a, &q);
        h.run(&a, 10);
        let expected = -0.5 * q.iter().map(|v| v * v).sum::<f64>();
        assert_eq!(h.dual_bound(), expected);
        // and the implied primal is the unconstrained minimizer
        assert_eq!(h.primal(), &[-1.0, 2.0, -0.5]);
    }

    /// Two variables, one constraint `u_1 - u_0 >= 0`, with `q = (-1, 1)` so the
    /// unconstrained minimizer `u = (1, -1)` violates it and the constraint is active.
    /// The exact answer pools both to their mean: `u_0 = u_1 = 0`, objective 0.
    #[test]
    fn active_constraint_pools_to_hand_derived_answer() {
        let edges = [(0usize, 1usize)];
        let coef = [(-1.0, 1.0)];
        let a = Constraints {
            edges: &edges,
            coef: &coef,
            n: 2,
        };
        let q = [-1.0, 1.0];

        let mut h = Hildreth::new(2, 1);
        h.seed_zero(&a, &q);
        h.run(&a, 200);

        let u = h.primal();
        assert!((u[0] - 0.0).abs() < 1e-12, "u = {u:?}");
        assert!((u[1] - 0.0).abs() < 1e-12, "u = {u:?}");
        assert!((h.dual_bound() - 0.0).abs() < 1e-12);
        assert!((h.multipliers()[0] - 1.0).abs() < 1e-12);
    }

    /// The bound must never exceed the primal value of a feasible point -- that is weak
    /// duality, and the whole reason the certificate is trustworthy. Checked at every
    /// sweep count, including zero, on a chain where the data is strictly anti-monotone
    /// so every constraint binds and the answer is the global mean.
    #[test]
    fn bound_never_exceeds_a_feasible_primal() {
        let n = 6;
        let edges: Vec<(usize, usize)> = (0..n - 1).map(|i| (i, i + 1)).collect();
        let coef: Vec<(f64, f64)> = vec![(-1.0, 1.0); n - 1];
        let a = Constraints {
            edges: &edges,
            coef: &coef,
            n,
        };
        // Strictly decreasing targets, so the isotonic answer pools everything.
        let q: Vec<f64> = (0..n).map(|i| i as f64 - 2.5).collect();

        // Feasible point: the constant vector at the mean of `-q`.
        let mean = -q.iter().sum::<f64>() / n as f64;
        let feasible = vec![mean; n];
        let primal_value = objective(&feasible, &q);

        let mut previous = f64::NEG_INFINITY;
        for sweeps in [0u32, 1, 2, 5, 20, 100, 1000] {
            let mut h = Hildreth::new(n, edges.len());
            h.seed_zero(&a, &q);
            h.run(&a, sweeps);
            let bound = h.dual_bound();
            assert!(
                bound <= primal_value + 1e-12,
                "weak duality violated at {sweeps} sweeps: {bound} > {primal_value}"
            );
            assert!(
                bound >= previous - 1e-12,
                "bound must be monotone in sweeps"
            );
            previous = bound;
        }
        // And it converges to the primal value, since the constant vector is optimal here.
        assert!(
            (previous - primal_value).abs() < 1e-9,
            "{previous} vs {primal_value}"
        );
    }
}
