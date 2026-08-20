pub mod admm;
pub mod hazard_rate_order;

use crate::partial_order::algorithm::uncensored::admm::constraints::Constraints;
use crate::partial_order::algorithm::uncensored::admm::hildreth::{self, Hildreth};
use crate::partial_order::algorithm::uncensored::admm::kkt::KktFactor;
use crate::partial_order::algorithm::uncensored::admm::polish::Polisher;
use crate::partial_order::algorithm::uncensored::admm::{polish, solver};
use crate::partial_order::algorithm::uncensored::hazard_rate_order::update_constraint_coefficients;
use crate::partial_order::routines::derive_transitive_reduction;
use crate::partial_order::{
    AlgorithmOutput, Config, OrderingInfo, QualityIndicators, UncensoredContext,
};
use crate::progress::ProgressTracker;
use crate::routines::transpose;
use crate::structures::{Direction, Increasing};
use crate::total_order::tonic_regression_pre_sorted;
use std::iter::repeat_n;

#[must_use]
pub fn algorithm<D: Direction, const HRO: bool>(
    context: &UncensoredContext<f64, f64>,
    config: Config,
    progress: &dyn ProgressTracker,
) -> AlgorithmOutput {
    // Build comparable pairs and reduce to cover edges
    let constraint_edges =
        derive_transitive_reduction(&context.x, context.n_covariate(), context.dimension());

    progress.set_total(context.n_threshold());

    if context.n_threshold() == 1 {
        // will always be Some
        // Single threshold -> all one's
        return AlgorithmOutput {
            cdfs: vec![1.0; context.n_covariate()],
            ordering_info: OrderingInfo::from_edges(constraint_edges, context.n_covariate()),
            quality_indicators: QualityIndicators {
                precision: 0.0,
                convergence_fraction: 1.0,
            },
        };
    }
    if context.n_covariate() == 1 {
        // Single covariate -> a single empirical cdf
        debug_assert!(constraint_edges.is_empty());
        return AlgorithmOutput {
            cdfs: context
                .weight
                .iter()
                .map(|w| w / context.x_weight[0])
                .scan(0.0, |acc, share| {
                    *acc += share;
                    Some(*acc as f32)
                })
                .collect(),
            ordering_info: OrderingInfo::from_edges(constraint_edges, context.n_covariate()),
            quality_indicators: QualityIndicators {
                precision: 0.0,
                convergence_fraction: 1.0,
            },
        };
    }

    // TODO: Remove covariates that are unconstraint

    // Build the coefficients of A, one pair per cover edge. We do an antitonic
    // (non-increasing) regression, so we reverse the order.
    let mut coefficients = build_order_coefficients::<D::REVERSE, HRO>(&constraint_edges);

    // --- Energy coordinates -------------------------------------------------------------
    // OSQP's termination tolerances are absolute on the problem it is handed (eps_abs, plus
    // eps_rel times the residual's own scale). The raw per-threshold QP has P = diag(w) and
    // q proportional to the weights, so both the global weight scale and within-fit weight
    // ratios silently weaken the enforced optimality: a coordinate with a tiny weight
    // contributes almost nothing to the dual residual and can stop far from its optimum
    // while OSQP reports a perfect solve. The estimator itself is invariant to rescaling
    // all weights by a positive constant, so none of this is statistically meaningful.
    // Two measures restore weight-scale robustness:
    //  1. Global power-of-two normalization `weight_scale` placing the largest group
    //     weight in (1/2, 1]. Power-of-two scaling is exact in binary floating point, so
    //     the solved problem is bitwise identical under rescaling all weights by a power
    //     of two (and equal up to rounding for any other positive factor).
    //  2. The change of variables u_i = sqrt_weight[i] * x_i with
    //     sqrt_weight[i] = sqrt(weight_scale * w_i). In u the objective is
    //     0.5 u^T u + q_energy^T u with q_energy[i] = weight_scale * q[i] / sqrt_weight[i],
    //     i.e. P is the identity (perfectly conditioned) and the per-coordinate dual
    //     residual scales like sqrt(w_i) instead of w_i: a uniform absolute tolerance now
    //     enforces per-coordinate optimality up to weight ratios of ~(1/eps)^2 instead
    //     of ~1/eps.
    // Constraint rows are mapped into u by dividing each column's coefficient by
    // sqrt_weight (see `scale_constraints_to_energy` and `update_constraint_matrix`); an
    // equivalent positive row rescaling keeps the largest |coefficient| at 1 so row values
    // stay in [-1, 1] for x in [0, 1]^n and the fixed row upper bound of 1.0 below stays
    // redundant.
    let max_weight = context.x_weight.iter().copied().fold(0.0, f64::max);
    debug_assert!(max_weight > 0.0 && max_weight.is_finite());
    // Clamp keeps 2^-exponent finite/normal even for absurd weight magnitudes.
    let mut weight_scale = 2f64.powi(-(max_weight.log2().ceil() as i32).clamp(-1022, 1022));
    if max_weight * weight_scale > 1.0 {
        // log2 + ceil can land just below a power-of-two boundary; one halving fixes it.
        weight_scale *= 0.5;
    }
    let sqrt_weight: Vec<f64> = context
        .x_weight
        .iter()
        .map(|w| (w * weight_scale).sqrt())
        .collect();
    scale_coefficients_to_energy(&mut coefficients, &constraint_edges, &sqrt_weight);

    // Construct the linear cost vector `q` in the original x-coordinates. It is the negative
    // of the element-wise product of (grouped) weights and mean response. Together with a
    // diagonal P = diag(w), it represents each component of the loss like (with p_i share of
    // weight below current threshold):
    // w_i (x_i - p_i)^2 = w_i (x_i^2 - 2 x_i p_i + p_i^2)
    //                   = w_i x_i^2 - 2 w_i x_i p_i + w_i p_i^2
    //                   ~ w_i x_i^2 - 2 w_i x_i p_i
    //                   ~ 0.5 w_i x_i^2 - w_i x_i p_i
    //                   = 0.5 x_i w_i x_i - w_i p_i x_i
    //                   = 0.5 x^T diag(w) x + q^T x
    // So to know the loss in the original space, double solver loss output and add sum(w_i p_i^2).
    //
    // In the above, 0 <= p_i = sum_j(w_ij / w_i if y_ij <= z) <= 1 and 0 < w_i = sum_j w_ij is the
    // total weight of the covariate, so we can equivalently use q_i = -sum_j(w_ij if y_ij <= z).
    //
    // The QP is solved in energy coordinates (see above): `q` keeps accumulating raw weights
    // across thresholds and `q_energy` is the mapped copy handed to the solver before each
    // solve. The accumulation stays in f64 regardless of the solver's working precision:
    // it is a running sum over every observation, so narrowing it would drift.
    //
    // Under HRO, we work in the reverse direction because we model the survival function.
    let mut q: Vec<f64> = if HRO {
        context.x_weight.iter().map(|w| -w).collect()
    } else {
        vec![0.0; context.n_covariate()]
    };
    // q_energy[i] = q[i] * weight_scale / sqrt_weight[i]; buffers reused across thresholds.
    let q_to_energy: Vec<f64> = sqrt_weight.iter().map(|sw| weight_scale / sw).collect();
    let mut q_energy = vec![0.0; context.n_covariate()];
    fn map_into(dst: &mut [f64], src: &[f64], factor: &[f64]) {
        for ((d, &s), &f) in dst.iter_mut().zip(src).zip(factor) {
            *d = s * f;
        }
    }

    let settings = solver::Settings {
        verbose: config.solver_settings.verbose,
        eps_abs: config.solver_settings.eps_abs,
        eps_rel: config.solver_settings.eps_rel,
        max_iter: config.solver_settings.max_iter,
    };
    // Only `n_threshold - 1` solves happen: at the last threshold every CDF is 1 (every
    // survival 0), which is filled in analytically below.
    let n_solve = context.n_threshold() - 1;
    let data_bounds = threshold_data_bounds(&context.y);
    debug_assert_eq!(data_bounds.len(), context.n_threshold() + 1);

    // We collect these results each iteration
    let mut cdfs = Vec::with_capacity(context.n_threshold() * context.n_covariate());
    let stats = if HRO {
        // The hazard-rate constraints are rebuilt from the previous threshold's fit, so the
        // sequence cannot be cut anywhere: threshold `t` is not even a well-posed problem
        // until `t - 1` has been solved.
        let mut solver_state = SolverState::new(context.n_covariate(), &constraint_edges);
        // The warm start is x = 1, the survival before any threshold, i.e. u = sqrt_weight.
        solver_state.workspace.x.copy_from_slice(&sqrt_weight);
        // The previous threshold's solution in x-coordinates, read as the survival S.
        let mut survival = vec![0.0; context.n_covariate()];
        for threshold in 0..n_solve {
            for observation in data_bounds[threshold]..data_bounds[threshold + 1] {
                q[context.x_indices[observation]] += context.weight[observation];
            }
            if threshold > 0 {
                // A's values change with S, so its factorization has to be refreshed. The
                // pattern does not change, so the symbolic phase is untouched.
                update_constraint_coefficients::<D>(
                    &mut coefficients,
                    &constraint_edges,
                    &survival,
                    &sqrt_weight,
                );
                solver_state.workspace.invalidate_factor();
            }
            map_into(&mut q_energy, &q, &q_to_energy);
            let constraints = Constraints {
                edges: &constraint_edges,
                coef: &coefficients,
                n: context.n_covariate(),
            };
            solver_state.solve_threshold(
                &constraints,
                &q_energy,
                &settings,
                &sqrt_weight,
                threshold,
                &mut survival,
            );
            cdfs.extend_from_slice(&survival);
            progress.increment();
        }
        solver_state.stats
    } else {
        // Under stochastic dominance `A` is the same matrix at every threshold and `q` is a
        // prefix sum over the response-sorted data, so only `q` moves from one threshold
        // to the next and the matrix handed to the solver never changes. Consecutive
        // problems are therefore near-identical, which is exactly what the warm start
        // carried in the solver state exploits: the factorization built for the first
        // solve serves the whole sequence.
        let mut solver_state = SolverState::new(context.n_covariate(), &constraint_edges);
        let mut fitted = vec![0.0; context.n_covariate()];
        let constraints = Constraints {
            edges: &constraint_edges,
            coef: &coefficients,
            n: context.n_covariate(),
        };
        for threshold in 0..n_solve {
            for observation in data_bounds[threshold]..data_bounds[threshold + 1] {
                q[context.x_indices[observation]] -= context.weight[observation];
            }
            map_into(&mut q_energy, &q, &q_to_energy);
            solver_state.solve_threshold(
                &constraints,
                &q_energy,
                &settings,
                &sqrt_weight,
                threshold,
                &mut fitted,
            );
            cdfs.extend_from_slice(&fitted);
            progress.increment();
        }
        solver_state.stats
    };

    // Finish the last trivial threshold
    cdfs.extend(repeat_n(if HRO { 0.0 } else { 1.0 }, context.n_covariate()));
    progress.increment();

    // Convert to covariate-major (a sequence of cdfs, one for each covariate)
    // TODO: This transpose can be avoided by computing the isotonic regressions first, then
    //  overwriting all at once (but this would require more memory for the regressions)
    transpose(&mut cdfs, context.n_threshold(), context.n_covariate());

    // From survival curve to CDF
    if HRO {
        for s in &mut cdfs {
            *s = 1.0 - *s;
        }
    }

    // Diagnostics before PAVA: precision = abs(min negative step across thresholds)
    let precision = cdfs
        .chunks_exact(context.n_threshold())
        .flat_map(|cdf| cdf.windows(2).map(|w| w[1] - w[0]))
        .reduce(f64::min) // min
        .map(|diff| if diff < 0.0 { -diff } else { 0.0 })
        .unwrap();

    // Apply PAVA along thresholds for each row, then append 1.0
    for cdf in cdfs.chunks_mut(context.n_threshold()) {
        let increasing =
            tonic_regression_pre_sorted::<Increasing, _, _>(cdf.iter().map(|&v| (v, 1.0)));
        for (initial, cleaned) in cdf.iter_mut().zip(increasing) {
            *initial = cleaned;
        }
    }

    if config.solver_settings.verbose {
        let SolveStats {
            unconverged,
            polished,
            iterations,
            refactorizations,
            rho_low,
            rho_high,
            f32_solves,
        } = stats;
        eprintln!(
            "solver: {n_solve} solves, {unconverged} unconverged, \
             {polished} polished, {iterations} iterations, \
             {refactorizations} refactorizations, rho in [{rho_low:.3e}, {rho_high:.3e}], \
             {f32_solves} solved with an f32 factor"
        );
    }
    // The last threshold is filled in analytically above, so it is not a solve and must not
    // dilute the fraction.
    let convergence_fraction = 1.0 - stats.unconverged as f64 / n_solve as f64;

    // The solver and the warm-started PAVA above run in f64; narrow once at the algorithm boundary
    // to match `AlgorithmOutput::cdfs`.
    AlgorithmOutput {
        cdfs: cdfs.into_iter().map(|v| v as f32).collect(),
        ordering_info: OrderingInfo::from_edges(constraint_edges, context.n_covariate()),
        quality_indicators: QualityIndicators {
            precision,
            convergence_fraction,
        },
    }
}

/// Where each threshold's observations sit in the response-sorted data: threshold `t` owns
/// `bounds[t]..bounds[t + 1]`.
///
/// The threshold sequence is the runs of equal responses, so one pass over `y` recovers it.
/// Having the boundaries up front keeps the accumulation of `q` a plain indexed walk
/// instead of a scan that has to re-detect where each threshold begins.
fn threshold_data_bounds(y: &[f64]) -> Vec<usize> {
    let mut bounds = vec![0usize];
    let mut observation = 0;
    while observation < y.len() {
        let value = y[observation];
        while observation < y.len() && y[observation] == value {
            observation += 1;
        }
        bounds.push(observation);
    }
    bounds
}

/// What the sequence of solves did, for diagnostics and the caller's convergence
/// accounting.
#[derive(Clone, Copy)]
struct SolveStats {
    unconverged: usize,
    polished: usize,
    iterations: u64,
    refactorizations: u64,
    rho_low: f64,
    rho_high: f64,
    f32_solves: usize,
}

impl SolveStats {
    fn new() -> Self {
        Self {
            unconverged: 0,
            polished: 0,
            iterations: 0,
            refactorizations: 0,
            rho_low: f64::INFINITY,
            rho_high: 0.0,
            f32_solves: 0,
        }
    }
}

/// Everything the threshold sequence solves with.
///
/// Sized once per fit and reused across thresholds. The pattern of the ADMM
/// system is the cover graph plus a diagonal and never changes, so the symbolic
/// factorization in `KktFactor::new` is the only one a block performs; its later thresholds
/// refresh values only. `workspace` carries the iterate and the adapted `rho` from one
/// threshold to the next, which is the warm start.
///
/// Bundled rather than kept as four locals so the threshold loop reads as the sequence of
/// solves it is, with the reused state named once.
struct SolverState {
    workspace: solver::Workspace,
    factor: KktFactor,
    polisher: Polisher,
    certifier: Hildreth,
    stats: SolveStats,
}

impl SolverState {
    fn new(n_covariate: usize, edges: &[(usize, usize)]) -> Self {
        Self {
            workspace: solver::Workspace::new(n_covariate, edges.len()),
            factor: KktFactor::new(n_covariate, edges),
            polisher: Polisher::new(n_covariate, edges.len()),
            certifier: Hildreth::new(n_covariate, edges.len()),
            stats: SolveStats::new(),
        }
    }

    /// Solve one threshold, writing the fit -- mapped back out of energy coordinates and
    /// clamped to `[0, 1]` -- into `out`.
    fn solve_threshold(
        &mut self,
        constraints: &Constraints<'_>,
        q_energy: &[f64],
        settings: &solver::Settings,
        sqrt_weight: &[f64],
        threshold: usize,
        out: &mut [f64],
    ) {
        let outcome = solver::solve(
            constraints,
            q_energy,
            settings,
            &mut self.workspace,
            &mut self.factor,
        );
        self.stats.iterations += u64::from(outcome.iterations);
        self.stats.refactorizations += u64::from(outcome.refactorizations);
        self.stats.rho_low = self.stats.rho_low.min(outcome.rho);
        self.stats.rho_high = self.stats.rho_high.max(outcome.rho);
        if outcome.used_f32 {
            self.stats.f32_solves += 1;
        }
        if outcome.status != solver::Status::Solved {
            self.stats.unconverged += 1;
            // Stalled. Dual coordinate ascent needs no factorization and is monotone in
            // the dual objective, so it can still make progress where the splitting has
            // stopped, and the point it implies is taken only if it is feasible and
            // strictly better -- the same gate the polish uses.
            self.certifier
                .seed_from_admm_dual(constraints, q_energy, &self.workspace.y);
            self.certifier.run(constraints, hildreth::STALL_SWEEPS);
            self.certifier
                .take_if_better(constraints, q_energy, &mut self.workspace.x);
            if settings.verbose {
                eprintln!(
                    "  threshold {threshold}: not converged (r_prim {:.2e}, r_dual {:.2e}), \
                     certified gap {:.2e}",
                    outcome.primal_residual,
                    outcome.dual_residual,
                    self.certifier.certified_gap(q_energy, &self.workspace.x),
                );
            }
        }

        // Polish to the exact minimizer over the identified active set. The gate inside
        // only accepts a feasible point with a strictly lower objective, so this can move
        // the answer towards the optimum but never away from it; the solver's own iterate
        // stands whenever the active set was not identified cleanly.
        let report = polish::polish(
            constraints,
            q_energy,
            &self.workspace.clipped,
            &self.workspace.x,
            &mut self.polisher,
        );
        if report.accepted {
            self.stats.polished += 1;
        }
        let primal = if report.accepted {
            &self.polisher.u
        } else {
            &self.workspace.x
        };

        // Map the primal solution back from energy coordinates: x_i = u_i / sqrt_weight[i].
        for ((value, &u), &sw) in out.iter_mut().zip(primal).zip(sqrt_weight) {
            *value = (u / sw).clamp(0.0, 1.0);
        }
    }
}

/// Build the coefficients of A from the cover edges (x_i <= x_j).
///
/// Row `r` encodes edge `constraints[r] = (i, j)` and carries exactly two nonzeros, -1 at
/// `i` and +1 at `j`, for the constraint `a^T x = x_j - x_i >= 0`. If HRO, the signs get
/// flipped, because then we're fitting a survival curve and not a CDF.
///
/// The returned pair is `(coefficient at i, coefficient at j)`; the edge list itself
/// supplies the column indices, so no sparse index structure is needed. The ±1
/// coefficients are in the original x-coordinates; `algorithm` rescales them into energy
/// coordinates with `scale_coefficients_to_energy` before the first solve.
pub fn build_order_coefficients<D: Direction, const HRO: bool>(
    constraints: &[(usize, usize)],
) -> Vec<(f64, f64)> {
    let coefficients = if D::IS_INCREASING != HRO {
        (-1.0, 1.0)
    } else {
        (1.0, -1.0)
    };
    vec![coefficients; constraints.len()]
}

/// Rescale the ±1 coefficients built by `build_order_coefficients` into energy coordinates
/// `u_i = sqrt_weight[i] * x_i`.
///
/// Expressing the row's constraint in `u` divides the column-`c` coefficient by
/// `sqrt_weight[c]`; multiplying the whole row by `min(sqrt_weight[i], sqrt_weight[j])` (a
/// positive row scaling, hence an equivalent constraint) then brings the largest
/// |coefficient| back to exactly 1. Keeping the coefficients bounded by 1 keeps every row
/// of the ADMM system on a common scale, so the solver's absolute residual tolerance means
/// the same thing on every constraint.
///
/// This is exactly the S = 1 case of the HRO rule in `update_constraint_coefficients`,
/// applied once up front: the SD coefficients are constant across thresholds, and the HRO
/// ones start from the plain order before the first per-threshold ratio update.
fn scale_coefficients_to_energy(
    coefficients: &mut [(f64, f64)],
    constraints: &[(usize, usize)],
    sqrt_weight: &[f64],
) {
    debug_assert_eq!(coefficients.len(), constraints.len());
    for (coefficient, &(i, j)) in coefficients.iter_mut().zip(constraints) {
        let scale = sqrt_weight[i].min(sqrt_weight[j]);
        coefficient.0 *= scale / sqrt_weight[i];
        coefficient.1 *= scale / sqrt_weight[j];
    }
}

#[cfg(test)]
mod test {
    use crate::IsotonicDistributionalRegressionFit;
    use crate::partial_order::structures::Fit;
    use crate::partial_order::{Config, CovariateGroups, Csr, OrderingInfo, SolverSettings};
    use crate::structures::StochasticOrder;

    /// Per-step survival ratios are nondecreasing along the covariate order. Same
    /// data and expectation as the total-order HRO kernel's doctest; the t=7 values
    /// are the KKT solution with both ratio constraints active.
    #[test]
    fn hazard_rate_ratio_constraint_orientation() {
        let fit: Fit<f64, f64> = Fit::fit::<f64>(
            &[1.0, 2.0, 3.0],
            &[8.0, 6.0, 7.0],
            None,
            None,
            CovariateGroups::empty(1),
            StochasticOrder::HazardRateOrder,
            false,
            Config {
                solver_settings: SolverSettings::default()
                    .verbose(false)
                    .eps_abs(1e-8)
                    .eps_rel(1e-8)
                    .max_iter(20_000),
            },
            &crate::NoProgress,
        )
        .unwrap();

        assert_eq!(fit.thresholds(), &[6.0, 7.0, 8.0]);
        #[rustfmt::skip]
        let expected = [
            0.5, 5.0 / 6.0, 1.0, // x = 1
            0.5, 5.0 / 6.0, 1.0, // x = 2
            0.0, 2.0 / 3.0, 1.0, // x = 3
        ];
        for (got, want) in fit.cdfs.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-3, "{:?}", fit.cdfs);
        }
    }

    fn fit_1d(
        x: &[f64],
        y: &[f64],
        w: &[f64],
        order: StochasticOrder,
        decreasing: bool,
    ) -> Fit<f64, f64> {
        Fit::fit(
            x,
            y,
            None,
            Some(w),
            CovariateGroups::empty(1),
            order,
            decreasing,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap()
    }

    /// HRO ratio constraints are oriented by the statically known edge direction. A
    /// pooled block makes the previous fitted survivals tied up to solver round-off,
    /// which must not decide the constraint direction: here t=0 pools all survivals
    /// to 6/7, and at t=1 the ratio constraints reduce to plain isotonicity, pooling
    /// to 4/7 — the same [1/7, 3/7, 1] the total-order kernel fits for every
    /// covariate.
    #[test]
    fn hazard_rate_pooled_block_keeps_orientation() {
        let fit = fit_1d(
            &[0.0, 1.0, 2.0, 2.0],
            &[1.5, 1.5, 0.0, 1.0],
            &[1.0, 3.0, 1.0, 2.0],
            StochasticOrder::HazardRateOrder,
            false,
        );
        assert_eq!(fit.thresholds(), &[0.0, 1.0, 1.5]);
        let expected = [1.0 / 7.0, 3.0 / 7.0, 1.0];
        for xv in [0.0, 1.0, 2.0] {
            let cdf: Vec<f32> = fit.cdf(&[xv]).collect();
            for (got, want) in cdf.iter().zip(expected.iter()) {
                assert!((got - want).abs() < 1e-3, "x = {xv}: {cdf:?}");
            }
        }
    }

    /// Decreasing direction with exactly tied fitted survivals: the static edge
    /// orientation keeps the empirical CDFs (which already satisfy the decreasing
    /// hazard order) as the fit.
    #[test]
    fn hazard_rate_decreasing_tied_survivals_keep_orientation() {
        let fit = fit_1d(
            &[1.0, 1.0, 2.0, 2.0],
            &[1.0, 3.0, 1.0, 2.0],
            &[1.0; 4],
            StochasticOrder::HazardRateOrder,
            true,
        );
        assert_eq!(fit.thresholds(), &[1.0, 2.0, 3.0]);
        let expected_x1 = [0.5, 0.5, 1.0];
        let expected_x2 = [0.5, 1.0, 1.0];
        for (xv, expected) in [(1.0, expected_x1), (2.0, expected_x2)] {
            let cdf: Vec<f32> = fit.cdf(&[xv]).collect();
            for (got, want) in cdf.iter().zip(expected.iter()) {
                assert!((got - want).abs() < 1e-3, "x = {xv}: {cdf:?}");
            }
        }
    }

    /// The QP is solved on power-of-two-normalized weights, so rescaling all weights
    /// by a power of two solves the bitwise-identical problem: the fits are equal.
    #[test]
    fn sd_fit_invariant_under_weight_rescaling() {
        let x = [1.0, 2.0, 3.0];
        let y = [2.0, 1.0, 3.0];
        let unit = fit_1d(
            &x,
            &y,
            &[1.0; 3],
            StochasticOrder::StochasticDominance,
            false,
        );
        let scaled = fit_1d(
            &x,
            &y,
            &[2f64.powi(-14); 3],
            StochasticOrder::StochasticDominance,
            false,
        );
        assert_eq!(unit.cdfs, scaled.cdfs);
    }

    /// Weight ratios within one fit must not degrade the solve: the QP runs in energy
    /// coordinates (u = sqrt(w)·x), so per-coordinate optimality holds for the tiny
    /// weight too. At z = 2 the empirical values [1, 1, 0] are already feasible, so
    /// the fitted row is the data.
    #[test]
    fn sd_fit_accurate_under_weight_ratio() {
        let fit = fit_1d(
            &[1.0, 2.0, 3.0],
            &[2.0, 1.0, 3.0],
            &[1024.0, 1.0, 0.0009765625],
            StochasticOrder::StochasticDominance,
            false,
        );
        assert_eq!(fit.thresholds(), &[1.0, 2.0, 3.0]);
        let x3: Vec<f32> = fit.cdf(&[3.0]).collect();
        let expected = [0.0, 0.0, 1.0];
        for (got, want) in x3.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-3, "x = 3: {x3:?}");
        }
        let f_x1_z1 = fit.cdf_at(&[1.0], 1.0);
        assert!(
            (f_x1_z1 - 1.0 / 1025.0).abs() < 1e-3,
            "F(x=1, z=1) = {f_x1_z1}, want 1/1025"
        );
    }

    #[test]
    fn multivariate_sd_case() {
        // 6 observations that collapse to 3 unique rows under SD preparation
        let covariates = [
            0.4, 0.6, 0.6, 0.4, // -> both become [0.6, 0.4]
            0.2, 0.3, 0.3, 0.2, // -> both become [0.3, 0.2]
            0.8, 0.9, 0.9, 0.8, // -> both become [0.9, 0.8]
        ];
        let responses = vec![0.2, 0.6, 0.1, 0.2, 0.7, 0.9];
        let weights = vec![1.0; responses.len()];
        let groups = CovariateGroups::parse([("sd", [0, 1])], 2).unwrap();

        let config = Config {
            solver_settings: SolverSettings::default()
                .verbose(false)
                .eps_abs(1e-8)
                .eps_rel(1e-8)
                .max_iter(10_000),
        };

        for so in [
            StochasticOrder::StochasticDominance,
            StochasticOrder::HazardRateOrder,
        ] {
            let fit = Fit::fit(
                &covariates,
                &responses,
                None,
                Some(&weights),
                groups.clone(),
                so,
                false,
                config.clone(),
                &crate::NoProgress,
            )
            .unwrap();
            fit.assert_consistent();

            assert!(fit.increasing);
            assert_eq!(fit.covariates.len(), 3 * 2);
            assert_eq!(fit.thresholds, vec![0.1, 0.2, 0.6, 0.7, 0.9]);
            if so == StochasticOrder::StochasticDominance {
                let expected = [1.0 / 6.0, 0.5, 2.0 / 3.0, 5.0 / 6.0, 1.0];
                for (&result, expected) in fit.global_cdf.iter().zip(expected.into_iter()) {
                    assert!((result - expected).abs() < 1e-6);
                }
                #[rustfmt::skip]
                let expected = [
                    0.5, 1.0, 1.0, 1.0, 1.0,
                    0.0, 0.5, 1.0, 1.0, 1.0,
                    0.0, 0.0, 0.0, 0.5, 1.0,
                ];
                for (&result, &expected) in fit.cdfs.iter().zip(expected.iter()) {
                    assert!((result - expected).abs() < 1e-6);
                }
            }
        }
    }

    #[test]
    fn multivariate_mixed_case() {
        // 6 observations -> 3 unique rows after SD (G1), ICX (G2), and COMP (G3) preparation.
        let covariates = vec![
            // a1, a2 | b1, b2 | c
            0.2, 0.1, 0.3, 0.2, 0.2, // R1a
            0.1, 0.2, 0.2, 0.3, 0.2, // R1b (permutes within G1 and G2)
            0.6, 0.4, 0.5, 0.3, 0.5, // R2a
            0.4, 0.6, 0.3, 0.5, 0.5, // R2b
            0.9, 0.8, 0.8, 0.7, 0.9, // R3a
            0.8, 0.9, 0.7, 0.8, 0.9, // R3b
        ];
        let responses = vec![0.1, 0.2, 0.2, 0.7, 0.6, 0.9];
        let weights = vec![1.0; responses.len()];

        // Groups: G1(a1,a2)=0 with SD, G2(b1,b2)=1 with ICX, G3(c)=2 with COMP
        let groups = [("sd", vec![0, 1]), ("icx", vec![2, 3]), ("comp", vec![4])];
        let groups = CovariateGroups::parse(groups, 5).unwrap();

        let config = Config {
            solver_settings: SolverSettings::default()
                .verbose(false)
                .eps_abs(1e-8)
                .eps_rel(1e-8)
                .max_iter(10_000),
        };

        let fit = Fit::fit(
            &covariates,
            &responses,
            None,
            Some(&weights),
            groups,
            StochasticOrder::StochasticDominance,
            false,
            config,
            &crate::NoProgress,
        )
        .unwrap();
        fit.assert_consistent();

        // This test assumes a specific sorting of the covariate values
        assert_eq!(
            fit.ordering_info,
            OrderingInfo {
                n: 3,
                smaller: Csr::from_sorted(&[(0, 1), (1, 2)], |(i, j)| (j, i), 3),
                larger: Csr::from_sorted(&[(0, 1), (1, 2)], |(i, j)| (i, j), 3),
                min: vec![0],
                max: vec![2],
            },
        );

        for (&result, expected) in fit
            .global_cdf
            .iter()
            .zip([1.0 / 6.0, 0.5, 2.0 / 3.0, 5.0 / 6.0, 1.0].into_iter())
        {
            assert!((result - expected).abs() < 1e-6);
        }

        assert_eq!(fit.thresholds, vec![0.1, 0.2, 0.6, 0.7, 0.9]);
        assert!(
            fit.quality_indicators.precision < 1e-6
                && fit.quality_indicators.convergence_fraction == 1.0
        );

        #[rustfmt::skip]
        let expected = [
            0.5, 1.0, 1.0, 1.0, 1.0, // R1
            0.0, 0.5, 0.5, 1.0, 1.0, // R2
            0.0, 0.0, 0.5, 0.5, 1.0, // R3
        ];
        for (&realized, &expected) in fit.cdfs.iter().zip(expected.iter()) {
            assert!((realized - expected).abs() < 1e-5);
        }

        let result = covariates
            .chunks_exact(5)
            // Filter every second row that gets merged
            .enumerate()
            .filter(|(i, _)| i % 2 == 0)
            .flat_map(|(_, v)| fit.cdf(v));
        for (realized, expected) in result.zip(expected) {
            assert!((realized - expected).abs() < 1e-5);
        }
    }

    // Build a multi-group, multi-threshold dataset exercising:
    // - duplicates created by permutations in the SD group (rows collapse)
    // - permutations in the ICX group that become identical after cumsum
    // - a comp group with mixed comparabilities
    // - multiple thresholds {0, 0.3, 0.7, 1}
    // - varied weights including a very small positive weight
    #[test]
    fn sd_multivariate_edge_cases() {
        // Columns: [SD0, SD1, COMP0, ICX0, ICX1, COMP1]
        let sd_pairs = [
            [0.2, 0.6], // duplicate with permutation below
            [0.6, 0.2],
            [0.9, 0.8], // duplicate with permutation below
            [0.8, 0.9],
            [0.5, 0.5],
            [0.4, 0.1],
        ];
        let comp_pairs = [[0.1, 0.9], [0.7, 0.3], [0.4, 0.5]];
        let icx_pairs = [
            [0.8, 0.3], // permutation below -> same after ICX cumsum
            [0.3, 0.8],
            [0.9, 0.1], // permutation below -> same after ICX cumsum
            [0.1, 0.9],
        ];

        let mut covariates = Vec::new();
        for (t, sp) in sd_pairs.iter().enumerate() {
            let cp = comp_pairs[t % comp_pairs.len()];
            let ip = icx_pairs[t % icx_pairs.len()];
            covariates.extend([sp[0], sp[1], cp[0], ip[0], ip[1], cp[1]]);
            // add a "paired" line that will collapse under SD/ICX preparation
            if t % 2 == 0 {
                covariates.extend([sp[1], sp[0], cp[0], ip[1], ip[0], cp[1]]);
            }
        }
        let covariate_dimension = 6;
        let n = covariates.len() / covariate_dimension;

        // responses with 4 unique thresholds
        let levels = [0.0, 0.3, 0.7, 1.0];
        let responses = (0..n)
            .map(|i| levels[i % levels.len()])
            .collect::<Vec<f64>>();

        // positive weights; include variety and a very small weight
        let mut weights = vec![1.0; n];
        if n >= 3 {
            weights[2] = 2.0;
        }
        if n >= 5 {
            weights[4] = 0.5;
        }
        weights[n - 1] = 1e-4;

        let groups = [("SD", [0, 1]), ("CoMp", [2, 5]), ("icx", [3, 4])];
        let groups = CovariateGroups::parse(groups, 6).unwrap();

        // more relaxed tolerances and many iterations so we (likely) solve to optimality
        let fit = Fit::fit(
            &covariates,
            &responses,
            None,
            Some(&weights),
            groups,
            StochasticOrder::StochasticDominance,
            false,
            Config {
                solver_settings: SolverSettings::default()
                    .verbose(false)
                    .eps_abs(1e-6)
                    .eps_rel(1e-6)
                    .max_iter(20_000),
            },
            &crate::NoProgress,
        )
        .unwrap();
        fit.assert_consistent();

        #[rustfmt::skip]
        let expected = [
            1.0, 1.0, 1.0, 1.0,
            0.0, 0.0, 0.5, 1.0,
            0.5, 1.0, 1.0, 1.0,
            0.0, 0.0, 1.0, 1.0,
            0.0, 1.0, 1.0, 1.0,
            1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0,
        ];
        for (&realized, &expected) in fit.cdfs.iter().zip(expected.iter()) {
            assert!((realized - expected).abs() < 1e-5);
        }

        let predicted = fit.cdf(&[f64::INFINITY; 6]);
        let expected = [0.0, 0.0, 1.0 / 3.0, 1.0];
        for (realized, expected) in predicted.zip(expected) {
            assert!((realized - expected).abs() < 1e-5);
        }
        let predicted = fit.cdf(&[f64::NEG_INFINITY; 6]);
        let expected = [1.0, 1.0, 1.0, 1.0];
        for (realized, expected) in predicted.zip(expected) {
            assert!((realized - expected).abs() < 1e-5);
        }
    }
}
