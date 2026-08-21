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
use crate::structures::Direction;
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
    // q_energy[i] = q[i] * weight_scale / sqrt_weight[i]; gathered per band below.
    let q_to_energy: Vec<f64> = sqrt_weight.iter().map(|sw| weight_scale / sw).collect();

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

    // --- Exact pin windows ----------------------------------------------------------------
    // Isotonic projection preserves bounds, so a covariate's fitted value is *exactly* 0 or
    // 1 outside a data-determined threshold window (see `pin_windows`). Outside its window
    // a covariate contributes nothing to the QP, and its constraint rows are vacuous, so
    // each run of thresholds is solved on the subgraph of in-window covariates only -- the
    // band where the CDF actually transitions. On response-correlated data that band is a
    // small fraction of the graph; on adversarial data it is everything and the reduction
    // degenerates to the full problem at no extra cost beyond the window sweep.
    //
    // Just as important as the size reduction: the pinned regions are exactly where the
    // fitted surface is flat, i.e. where constraint rows are degenerate (margin zero) and
    // the active set is ambiguous. Excluding them keeps the solver's active-set estimate --
    // and with it the exact-finish certificates -- sharp.
    let (miny, maxy) = covariate_response_extremes(context);
    let windows = pin_windows(
        &constraint_edges,
        context.n_covariate(),
        &miny,
        &maxy,
        D::IS_INCREASING,
    );

    let n_covariate = context.n_covariate();
    let n_edges = constraint_edges.len();
    // We collect these results each iteration
    let mut cdfs = Vec::with_capacity(context.n_threshold() * n_covariate);
    let mut stats = SolveStats::new();

    // Warm-start carriers across band rebuilds, in the full index space: the primal in
    // energy coordinates, the dual per cover edge, and the active-set record per cover
    // edge. A covariate entering the band comes from its pinned state (value zero, dual
    // zero on its edges), which is exactly the fresh entries' initialization.
    let mut full_x = if HRO {
        // The warm start is x = 1, the survival before any threshold, i.e. u = sqrt_weight.
        sqrt_weight.clone()
    } else {
        vec![0.0; n_covariate]
    };
    let mut full_y = vec![0.0; n_edges];
    let mut full_clip = vec![false; n_edges];
    // The previously *emitted* row in x-coordinates. It anchors the threshold-direction
    // clip below, and on the hazard-rate path it is also what the next threshold's ratio
    // constraints are rebuilt from -- which is why that path cannot cut the threshold
    // sequence: threshold `t` is not even a well-posed problem until `t - 1` has been
    // solved. Before any threshold the CDF is 0 everywhere (the survival 1).
    let mut previous_fit = vec![if HRO { 1.0 } else { 0.0 }; n_covariate];
    let mut row = vec![0.0; n_covariate];
    // Largest monotonicity violation the composition below had to absorb, reported as
    // `QualityIndicators::precision`.
    let mut precision = 0.0f64;

    let mut band: Option<(SubInstance, SolverState)> = None;
    let mut chunk_start = 0usize;
    while chunk_start < n_solve {
        let chunk_end = (chunk_start + CHUNK_THRESHOLDS).min(n_solve);
        let z_first = context.thresholds[chunk_start];
        let z_last = context.thresholds[chunk_end - 1];

        // In-window at some threshold of this chunk. The hazard-rate path has no lower
        // pin (survival = 1 has no closed-form certificate under ratio constraints), so
        // only the upper window bounds it.
        let active: Vec<usize> = (0..n_covariate)
            .filter(|&c| windows.one[c] > z_first && (HRO || windows.zero[c] <= z_last))
            .collect();

        // Rebuild the band only when the active set moved; while it is unchanged the
        // existing state -- symbolic factorization included -- continues seamlessly. On
        // data where nothing pins, this collapses to one build for the whole fit.
        let moved = match &band {
            Some((instance, _)) => instance.active != active,
            None => !active.is_empty(),
        };
        if moved {
            if let Some((instance, state)) = band.take() {
                instance.scatter_state(&state, &mut full_x, &mut full_y, &mut full_clip);
                stats.merge(state.stats);
            }
            if !active.is_empty() {
                let instance = SubInstance::build(
                    active,
                    n_covariate,
                    &constraint_edges,
                    &coefficients,
                    &sqrt_weight,
                    &q_to_energy,
                );
                let mut state = SolverState::new(instance.active.len(), &instance.edges);
                instance.gather_state(&mut state, &full_x, &full_y, &full_clip);
                band = Some((instance, state));
            }
        }

        for threshold in chunk_start..chunk_end {
            let z = context.thresholds[threshold];
            for observation in data_bounds[threshold]..data_bounds[threshold + 1] {
                let covariate = context.x_indices[observation];
                if HRO {
                    q[covariate] += context.weight[observation];
                } else {
                    q[covariate] -= context.weight[observation];
                }
            }

            if let Some((instance, state)) = band.as_mut() {
                if HRO && threshold > 0 {
                    // A's values change with S, so its factorization has to be refreshed.
                    // The pattern does not change, so the symbolic phase is untouched.
                    for (slot, &covariate) in instance.scratch.iter_mut().zip(&instance.active) {
                        *slot = previous_fit[covariate];
                    }
                    update_constraint_coefficients::<D>(
                        &mut instance.coefficients,
                        &instance.edges,
                        &instance.scratch,
                        &instance.sqrt_weight,
                    );
                    state.workspace.invalidate_factor();
                }
                for (slot, (&covariate, &scale)) in instance
                    .q_energy
                    .iter_mut()
                    .zip(instance.active.iter().zip(&instance.q_scale))
                {
                    *slot = q[covariate] * scale;
                }
                let constraints = Constraints {
                    edges: &instance.edges,
                    coef: &instance.coefficients,
                    n: instance.active.len(),
                };
                let mut fitted = std::mem::take(&mut instance.fitted);
                state.solve_threshold(
                    &constraints,
                    &instance.q_energy,
                    &settings,
                    &instance.sqrt_weight,
                    threshold,
                    &mut fitted,
                );
                instance.fitted = fitted;
            }

            // Assemble the output row: exact pins outside the window, solved values
            // inside. On the hazard-rate path the pinned zeros feed the next threshold's
            // ratio constraints as exact zeros, which is what routes those rows onto the
            // plain-order fallback.
            for (c, (value, (&one_c, &zero_c))) in row
                .iter_mut()
                .zip(windows.one.iter().zip(&windows.zero))
                .enumerate()
            {
                *value = if z >= one_c {
                    if HRO { 0.0 } else { 1.0 }
                } else if !HRO && z < zero_c {
                    0.0
                } else {
                    let (instance, _) = band
                        .as_ref()
                        .expect("an in-window covariate must be inside the band");
                    instance.fitted[instance.map[c]]
                };
            }

            // --- Inherent monotonicity, both directions -------------------------------
            // The exact optimum satisfies the covariate order and is monotone in the
            // threshold, but a solve only enforces either to tolerance. Both properties
            // are therefore established here, on the emitted values themselves, the same
            // way PAVA earns its guarantee: by construction.
            //
            //  1. One sweep over the band's cover edges in topological order pins every
            //     edge inequality exactly: the read endpoint is final before any edge
            //     leaving it is processed, and the written endpoint only moves further in
            //     the enforced direction. Transitivity then covers all comparable pairs,
            //     and edges outside the band hold already -- a pinned endpoint's
            //     inequality is vacuous for values in [0, 1], and the pin windows are
            //     closed in exactly the directions that make pins mutually consistent.
            //  2. Clipping against the previously emitted row makes the threshold
            //     direction exact, and cannot break step 1: the previous row satisfies
            //     the same edge inequalities by induction, and min/max are monotone in
            //     both arguments, so the clipped row inherits the order from its two
            //     ordered inputs.
            //
            // Both moves act only on solver noise -- the exact solution is a fixed point
            // -- and the largest absorbed violation is reported as the fit's `precision`.
            if let Some((instance, _)) = band.as_ref() {
                // The fit direction decides which way each cover edge points: the CDF is
                // antitone in the covariate order for an increasing fit (the survival
                // antitone for a decreasing one), so exactly one of `min` and `max`
                // repairs toward feasibility.
                if D::IS_INCREASING != HRO {
                    for &(source, target) in &instance.sweep_edges {
                        let step = row[target] - row[source];
                        if step > 0.0 {
                            if step > precision {
                                precision = step;
                            }
                            row[target] = row[source];
                        }
                    }
                } else {
                    for &(source, target) in &instance.sweep_edges {
                        let step = row[source] - row[target];
                        if step > 0.0 {
                            if step > precision {
                                precision = step;
                            }
                            row[target] = row[source];
                        }
                    }
                }
            }
            for (value, &previous) in row.iter_mut().zip(previous_fit.iter()) {
                // CDFs may not fall across thresholds; survivals may not rise.
                let step = if HRO {
                    *value - previous
                } else {
                    previous - *value
                };
                if step > 0.0 {
                    if step > precision {
                        precision = step;
                    }
                    *value = previous;
                }
            }
            previous_fit.copy_from_slice(&row);
            cdfs.extend_from_slice(&row);
            progress.increment();
        }
        chunk_start = chunk_end;
    }
    if let Some((_, state)) = band.take() {
        stats.merge(state.stats);
    }
    let stats = stats;

    // Finish the last trivial threshold
    cdfs.extend(repeat_n(if HRO { 0.0 } else { 1.0 }, context.n_covariate()));
    progress.increment();

    // Convert to covariate-major (a sequence of cdfs, one for each covariate)
    transpose(&mut cdfs, context.n_threshold(), context.n_covariate());

    // From survival curve to CDF. `1 - s` and the later f64 -> f32 narrowing are both
    // monotone maps, so the two monotonicity guarantees established row by row above
    // survive them exactly.
    if HRO {
        for s in &mut cdfs {
            *s = 1.0 - *s;
        }
    }

    if config.solver_settings.verbose {
        let SolveStats {
            unconverged,
            polished,
            exact_finishes,
            iterations,
            refactorizations,
            rho_low,
            rho_high,
            f32_solves,
        } = stats;
        eprintln!(
            "solver: {n_solve} solves, {unconverged} unconverged, \
             {exact_finishes} exact finishes, {polished} polished, \
             {iterations} iterations, \
             {refactorizations} refactorizations, rho in [{rho_low:.3e}, {rho_high:.3e}], \
             {f32_solves} solved with an f32 factor"
        );
    }
    // The last threshold is filled in analytically above, so it is not a solve and must not
    // dilute the fraction.
    let convergence_fraction = 1.0 - stats.unconverged as f64 / n_solve as f64;

    // The solver runs in f64; narrow once at the algorithm boundary to match
    // `AlgorithmOutput::cdfs`.
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

/// Thresholds sharing one band of in-window covariates.
///
/// The pin windows change with every threshold, but rebuilding the band -- and with it
/// the solver state and its symbolic factorization -- per threshold would forfeit the
/// warm start. A run of this many thresholds shares the union of its windows instead:
/// covariates whose window opens mid-run are simply solved as free variables while their
/// exact pin overrides the output, so the fit is unchanged and only the reduction is
/// coarser. Sixteen keeps the union close to the per-threshold band (a window spans many
/// thresholds on any realistic response distribution) while amortizing a rebuild to a
/// sixteenth of a solve each.
const CHUNK_THRESHOLDS: usize = 16;

/// Smallest and largest response observed at each covariate.
fn covariate_response_extremes(context: &UncensoredContext<f64, f64>) -> (Vec<f64>, Vec<f64>) {
    let n = context.n_covariate();
    let mut miny = vec![f64::INFINITY; n];
    let mut maxy = vec![f64::NEG_INFINITY; n];
    for (&covariate, &response) in context.x_indices.iter().zip(&context.y) {
        if response < miny[covariate] {
            miny[covariate] = response;
        }
        if response > maxy[covariate] {
            maxy[covariate] = response;
        }
    }
    (miny, maxy)
}

/// The exact pin thresholds of every covariate.
///
/// `one[c]` and `zero[c]` bound the window outside which covariate `c`'s fitted value is
/// known in closed form, threshold by threshold:
///
///  * for `z >= one[c]` the fitted CDF is exactly 1 (fitted survival exactly 0 under the
///    hazard-rate order);
///  * for `z < zero[c]` the fitted CDF is exactly 0 (stochastic dominance only; the
///    hazard-rate ratio constraints admit no such closed form on the survival-1 side).
///
/// The argument, for an increasing fit (CDF antitone in the covariate order, `x_i >= x_j`
/// on every cover edge `i` before `j`): let `L` be the lower closure of `c`. If every
/// observation of every covariate in `L` has response `<= z`, take the optimum `x*` and
/// raise it to 1 on `L`. The result is still feasible -- `L` is lower-closed, so every
/// edge into `L` comes from inside `L`, and edges leaving `L` only see their lower side
/// raised to the maximum -- and its objective cannot be worse, since each raised
/// coordinate reaches its data value exactly. Strict convexity then forces `x* = 1` on
/// `L`. Hence `one[c] = max` of the response maxima over the lower closure, and by the
/// mirrored argument `zero[c] = min` of the response minima over the upper closure; a
/// decreasing fit (CDF isotone) swaps the closures. Both sweeps propagate along the cover
/// edges, whose endpoints are topologically ordered (`i < j`), in `O(n + m log m)`.
///
/// The bounds are exact, not conservative: the projection onto the monotone cone maps
/// `[0, 1]`-valued data to `[0, 1]`-valued fits, which is what makes the raised point
/// feasible. That is also why the pinned values can be written into the output verbatim
/// -- they are *more* accurate than anything an iterative solve would produce.
struct PinWindows {
    /// Fit is exactly 1 (CDF) / 0 (survival) for thresholds at or above this.
    one: Vec<f64>,
    /// Fit is exactly 0 (CDF) for thresholds strictly below this.
    zero: Vec<f64>,
}

fn pin_windows(
    edges: &[(usize, usize)],
    n: usize,
    miny: &[f64],
    maxy: &[f64],
    increasing: bool,
) -> PinWindows {
    debug_assert_eq!(miny.len(), n);
    debug_assert_eq!(maxy.len(), n);
    let mut one = maxy.to_vec();
    let mut zero = miny.to_vec();
    // Propagation along a cover edge must see its source's final value, so the edges are
    // walked in topological order of the propagation direction; `i < j` on every edge
    // makes sorting by the receiving endpoint sufficient.
    let mut order: Vec<usize> = (0..edges.len()).collect();
    if increasing {
        // `one` accumulates over lower closures, `zero` over upper closures.
        order.sort_unstable_by_key(|&r| edges[r].1);
        for &r in &order {
            let (i, j) = edges[r];
            one[j] = one[j].max(one[i]);
        }
        order.sort_unstable_by_key(|&r| std::cmp::Reverse(edges[r].0));
        for &r in &order {
            let (i, j) = edges[r];
            zero[i] = zero[i].min(zero[j]);
        }
    } else {
        // Mirrored: `one` over upper closures, `zero` over lower closures.
        order.sort_unstable_by_key(|&r| std::cmp::Reverse(edges[r].0));
        for &r in &order {
            let (i, j) = edges[r];
            one[i] = one[i].max(one[j]);
        }
        order.sort_unstable_by_key(|&r| edges[r].1);
        for &r in &order {
            let (i, j) = edges[r];
            zero[j] = zero[j].min(zero[i]);
        }
    }
    PinWindows { one, zero }
}

/// One band's quadratic program: the full problem restricted to the in-window covariates,
/// in a compact index space.
///
/// Restriction is exact, not approximate: a pinned covariate's rows are vacuous (a pin at
/// 1 only relaxes `x_i >= x_j` rows below it, a pin at 0 only rows above, and the pinned
/// sets are closed in the respective direction, so no row ever constrains a free
/// covariate through a pinned one), and its objective term is constant. The reduced
/// optimum therefore *is* the full optimum on the band.
struct SubInstance {
    /// Band index -> covariate.
    active: Vec<usize>,
    /// Covariate -> band index, `usize::MAX` outside the band.
    map: Vec<usize>,
    /// Cover edges with both endpoints in the band, in band indices.
    edges: Vec<(usize, usize)>,
    /// The global row each band row came from, for the dual warm-start carriers.
    edge_ids: Vec<usize>,
    coefficients: Vec<(f64, f64)>,
    sqrt_weight: Vec<f64>,
    /// `q_to_energy` restricted to the band.
    q_scale: Vec<f64>,
    q_energy: Vec<f64>,
    fitted: Vec<f64>,
    scratch: Vec<f64>,
    /// The band's cover edges as full-space covariate pairs, ordered by the receiving
    /// endpoint so the monotonicity sweep reads only finalized values. Every cover edge
    /// with a pinned endpoint holds exactly without repair -- pins are closed in the
    /// consistent directions and dominate any value in [0, 1] -- so this list is also the
    /// complete set of edges the sweep ever has to visit.
    sweep_edges: Vec<(usize, usize)>,
}

impl SubInstance {
    fn build(
        active: Vec<usize>,
        n: usize,
        edges: &[(usize, usize)],
        coefficients: &[(f64, f64)],
        sqrt_weight: &[f64],
        q_to_energy: &[f64],
    ) -> Self {
        let mut map = vec![usize::MAX; n];
        for (k, &covariate) in active.iter().enumerate() {
            map[covariate] = k;
        }
        let mut sub_edges = Vec::new();
        let mut edge_ids = Vec::new();
        let mut sub_coefficients = Vec::new();
        let mut sweep_edges = Vec::new();
        for (r, &(i, j)) in edges.iter().enumerate() {
            if map[i] != usize::MAX && map[j] != usize::MAX {
                sub_edges.push((map[i], map[j]));
                edge_ids.push(r);
                sub_coefficients.push(coefficients[r]);
                sweep_edges.push((i, j));
            }
        }
        // Topological for the sweep: `i < j` on every cover edge, so ordering by the
        // receiving endpoint guarantees all edges into a node precede all edges out of it.
        sweep_edges.sort_unstable_by_key(|&(_, j)| j);
        let k = active.len();
        Self {
            sqrt_weight: active.iter().map(|&c| sqrt_weight[c]).collect(),
            q_scale: active.iter().map(|&c| q_to_energy[c]).collect(),
            q_energy: vec![0.0; k],
            fitted: vec![0.0; k],
            scratch: vec![0.0; k],
            active,
            map,
            edges: sub_edges,
            edge_ids,
            coefficients: sub_coefficients,
            sweep_edges,
        }
    }

    /// Load the warm start from the full-space carriers into a fresh state.
    ///
    /// The slack is rebuilt as `max(0, A x)` rather than carried: it is determined by the
    /// primal up to the projection, and the band's edge set just changed.
    fn gather_state(&self, state: &mut SolverState, x: &[f64], y: &[f64], clip: &[bool]) {
        let ws = &mut state.workspace;
        for (slot, &covariate) in ws.x.iter_mut().zip(&self.active) {
            *slot = x[covariate];
        }
        for (slot, &row) in ws.y.iter_mut().zip(&self.edge_ids) {
            *slot = y[row];
        }
        for (slot, &row) in ws.clipped.iter_mut().zip(&self.edge_ids) {
            *slot = clip[row];
        }
        let constraints = Constraints {
            edges: &self.edges,
            coef: &self.coefficients,
            n: self.active.len(),
        };
        constraints.mul(&ws.x, &mut ws.z);
        for z_r in ws.z.iter_mut() {
            *z_r = z_r.max(0.0);
        }
    }

    /// Store the state back into the full-space carriers for the next band to pick up.
    fn scatter_state(&self, state: &SolverState, x: &mut [f64], y: &mut [f64], clip: &mut [bool]) {
        let ws = &state.workspace;
        for (&value, &covariate) in ws.x.iter().zip(&self.active) {
            x[covariate] = value;
        }
        for (&value, &row) in ws.y.iter().zip(&self.edge_ids) {
            y[row] = value;
        }
        for (&value, &row) in ws.clipped.iter().zip(&self.edge_ids) {
            clip[row] = value;
        }
    }
}

/// What the sequence of solves did, for diagnostics and the caller's convergence
/// accounting.
#[derive(Clone, Copy)]
struct SolveStats {
    unconverged: usize,
    polished: usize,
    /// Solves that ended on a certified duality gap instead of the residual test.
    exact_finishes: usize,
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
            exact_finishes: 0,
            iterations: 0,
            refactorizations: 0,
            rho_low: f64::INFINITY,
            rho_high: 0.0,
            f32_solves: 0,
        }
    }

    /// Fold one band's stats into the fit total. Every field combines by sum, minimum or
    /// maximum, so the total does not depend on where the band boundaries fell.
    fn merge(&mut self, other: Self) {
        self.unconverged += other.unconverged;
        self.polished += other.polished;
        self.exact_finishes += other.exact_finishes;
        self.iterations += other.iterations;
        self.refactorizations += other.refactorizations;
        self.rho_low = self.rho_low.min(other.rho_low);
        self.rho_high = self.rho_high.max(other.rho_high);
        self.f32_solves += other.f32_solves;
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
            &mut self.polisher,
        );
        self.stats.iterations += u64::from(outcome.iterations);
        self.stats.refactorizations += u64::from(outcome.refactorizations);
        self.stats.rho_low = self.stats.rho_low.min(outcome.rho);
        self.stats.rho_high = self.stats.rho_high.max(outcome.rho);
        if outcome.used_f32 {
            self.stats.f32_solves += 1;
        }
        if outcome.finished_exact {
            self.stats.exact_finishes += 1;
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
        // stands whenever the active set was not identified cleanly. A solve that finished
        // through a certified polish already returned exactly this point.
        let primal = if outcome.finished_exact {
            &self.workspace.x
        } else {
            let report = polish::polish(
                constraints,
                q_energy,
                |r| self.workspace.clipped[r],
                &self.workspace.x,
                &mut self.polisher,
            );
            // Whatever the acceptance verdict, leave the *multiplier support* of the
            // polished point as the recorded active set, not the raw projection clips.
            // Near-degenerate rows -- margin at rounding level, dual near zero -- get
            // clipped spuriously and glue pools together, and the next threshold's
            // predictor would then have to undo every spurious merge one exchange at a
            // time. The peel prices the rows exactly, so its support is the clean
            // structural summary; a row it drops wrongly is re-added by the repair loop
            // on the next attempt.
            self.polisher.peel_multipliers(constraints, q_energy, None);
            for (clip_r, &lambda_r) in self
                .workspace
                .clipped
                .iter_mut()
                .zip(&self.polisher.multipliers)
            {
                *clip_r = lambda_r > 0.0;
            }
            if report.accepted {
                self.stats.polished += 1;
                &self.polisher.u
            } else {
                &self.workspace.x
            }
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

    /// Deterministic 64-bit LCG, so the random sweeps below are reproducible.
    struct Lcg(u64);

    impl Lcg {
        fn next_u64(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0
        }

        /// Uniform in [0, 1).
        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }

        fn below(&mut self, bound: usize) -> usize {
            (self.next_u64() >> 33) as usize % bound
        }
    }

    /// `pin_windows` against the definition: the closure extremes computed by brute-force
    /// dominance over random 2-D posets, in both fit directions.
    #[test]
    fn pin_windows_match_bruteforce_closures() {
        use crate::partial_order::routines::derive_transitive_reduction;

        for seed in 0..50u64 {
            let mut rng = Lcg(0x0008_17d0_5eed ^ seed);
            let levels = 2 + rng.below(4) as i64;
            let mut points: Vec<(i64, i64)> = (0..(2 + rng.below(12)))
                .map(|_| {
                    (
                        rng.below(levels as usize) as i64,
                        rng.below(levels as usize) as i64,
                    )
                })
                .collect();
            points.sort_unstable();
            points.dedup();
            let n = points.len();
            let flat: Vec<f64> = points
                .iter()
                .flat_map(|&(a, b)| [a as f64, b as f64])
                .collect();
            let edges = derive_transitive_reduction(&flat, n, 2);
            let miny: Vec<f64> = (0..n).map(|_| rng.unit()).collect();
            let maxy: Vec<f64> = miny.iter().map(|&low| low + rng.unit()).collect();

            let below = |a: (i64, i64), b: (i64, i64)| a.0 <= b.0 && a.1 <= b.1;
            for &increasing in &[true, false] {
                let windows = super::pin_windows(&edges, n, &miny, &maxy, increasing);
                for c in 0..n {
                    // The closure feeding `one` for an increasing fit is the lower one.
                    let mut expected_one = f64::NEG_INFINITY;
                    let mut expected_zero = f64::INFINITY;
                    for j in 0..n {
                        let in_lower = below(points[j], points[c]);
                        let in_upper = below(points[c], points[j]);
                        if (increasing && in_lower) || (!increasing && in_upper) {
                            expected_one = expected_one.max(maxy[j]);
                        }
                        if (increasing && in_upper) || (!increasing && in_lower) {
                            expected_zero = expected_zero.min(miny[j]);
                        }
                    }
                    assert_eq!(
                        windows.one[c], expected_one,
                        "seed {seed}, increasing {increasing}, covariate {c}"
                    );
                    assert_eq!(
                        windows.zero[c], expected_zero,
                        "seed {seed}, increasing {increasing}, covariate {c}"
                    );
                }
            }
        }
    }

    /// Outside its pin window a covariate's fitted CDF is written verbatim as 0 or 1 --
    /// bitwise, through the public fit path, not merely within tolerance.
    #[test]
    fn pinned_regions_are_bitwise_exact() {
        let n = 60;
        let mut rng = Lcg(0x91d_2024);
        let mut covariates = Vec::with_capacity(2 * n);
        let mut responses = Vec::with_capacity(n);
        for _ in 0..n {
            let a = rng.unit();
            let b = rng.unit();
            covariates.push(a);
            covariates.push(b);
            responses.push(a + b + 0.5 * (rng.unit() - 0.5));
        }
        let groups = CovariateGroups::parse([("sd", [0usize, 1])], 2).unwrap();
        let fit: Fit<f64, f64> = Fit::fit::<f64>(
            &covariates,
            &responses,
            None,
            None,
            groups,
            StochasticOrder::StochasticDominance,
            false,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();

        let n_threshold = fit.thresholds.len();
        let n_covariate = fit.covariates.len() / 2;
        // Recover each fit covariate's response from the SD-prepared row (sorted
        // descending within the group), then take closure extremes by brute force.
        let prepared: Vec<(f64, f64)> = (0..n)
            .map(|observation| {
                let a = covariates[2 * observation];
                let b = covariates[2 * observation + 1];
                (a.max(b), a.min(b))
            })
            .collect();
        let rows: Vec<(f64, f64)> = (0..n_covariate)
            .map(|c| (fit.covariates[2 * c], fit.covariates[2 * c + 1]))
            .collect();
        let response_of = |row: (f64, f64)| {
            responses
                .iter()
                .zip(&prepared)
                .find(|&(_, &p)| p == row)
                .map(|(&y, _)| y)
                .expect("every fit covariate row comes from an observation")
        };
        for (c, &row) in rows.iter().enumerate() {
            // Increasing fit: `one` over the lower closure, `zero` over the upper.
            let mut one = f64::NEG_INFINITY;
            let mut zero = f64::INFINITY;
            for &other in &rows {
                if other.0 <= row.0 && other.1 <= row.1 {
                    one = one.max(response_of(other));
                }
                if other.0 >= row.0 && other.1 >= row.1 {
                    zero = zero.min(response_of(other));
                }
            }
            for (t, &z) in fit.thresholds.iter().enumerate() {
                let value = fit.cdfs[c * n_threshold + t];
                if z >= one {
                    assert_eq!(value, 1.0f32, "covariate {c}, threshold {z}");
                } else if z < zero {
                    assert_eq!(value, 0.0f32, "covariate {c}, threshold {z}");
                }
            }
        }
    }

    /// A decreasing fit must be the increasing fit of the negated covariates -- the two
    /// runs exercise the mirrored pin-window sweeps and constraint orientations against
    /// each other.
    #[test]
    fn decreasing_mirrors_negated_increasing() {
        let n = 40;
        let mut rng = Lcg(0xdec_2024);
        let mut covariates = Vec::with_capacity(2 * n);
        let mut negated = Vec::with_capacity(2 * n);
        let mut responses = Vec::with_capacity(n);
        for _ in 0..n {
            let a = rng.unit();
            let b = rng.unit();
            covariates.extend([a, b]);
            negated.extend([-a, -b]);
            responses.push(a + b + 0.5 * (rng.unit() - 0.5));
        }
        let groups = CovariateGroups::parse([("sd", [0usize, 1])], 2).unwrap();
        let decreasing: Fit<f64, f64> = Fit::fit::<f64>(
            &covariates,
            &responses,
            None,
            None,
            groups.clone(),
            StochasticOrder::StochasticDominance,
            true,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();
        let increasing: Fit<f64, f64> = Fit::fit::<f64>(
            &negated,
            &responses,
            None,
            None,
            groups,
            StochasticOrder::StochasticDominance,
            false,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();

        assert_eq!(decreasing.thresholds, increasing.thresholds);
        for observation in 0..n {
            let point = &covariates[2 * observation..2 * observation + 2];
            let mirrored = [-point[0], -point[1]];
            let from_decreasing: Vec<f32> = decreasing.cdf(point).collect();
            let from_increasing: Vec<f32> = increasing.cdf(&mirrored).collect();
            for (got, want) in from_decreasing.iter().zip(&from_increasing) {
                assert!(
                    (got - want).abs() < 1e-4,
                    "observation {observation}: {from_decreasing:?} vs {from_increasing:?}"
                );
            }
        }
    }

    /// The emitted CDFs are exactly monotone in both directions -- bitwise, not within a
    /// tolerance -- for every combination of stochastic order and fit direction. Checked
    /// over *all* comparable covariate pairs, not just cover edges, since violations
    /// would otherwise accumulate along chains.
    #[test]
    fn fits_are_bitwise_monotone_in_both_directions() {
        let n = 60;
        let mut rng = Lcg(0x2020_0821);
        let mut covariates = Vec::with_capacity(2 * n);
        let mut responses = Vec::with_capacity(n);
        for _ in 0..n {
            let a = rng.unit();
            let b = rng.unit();
            covariates.extend([a, b]);
            responses.push(a + b + 0.5 * (rng.unit() - 0.5));
        }
        let groups = CovariateGroups::parse([("sd", [0usize, 1])], 2).unwrap();

        for (label, order) in [
            ("sd", StochasticOrder::StochasticDominance),
            ("hro", StochasticOrder::HazardRateOrder),
        ] {
            for decreasing in [false, true] {
                let fit: Fit<f64, f64> = Fit::fit::<f64>(
                    &covariates,
                    &responses,
                    None,
                    None,
                    groups.clone(),
                    order,
                    decreasing,
                    Config::default(),
                    &crate::NoProgress,
                )
                .unwrap();
                let n_threshold = fit.thresholds.len();
                let n_covariate = fit.covariates.len() / 2;

                for c in 0..n_covariate {
                    let cdf = &fit.cdfs[c * n_threshold..(c + 1) * n_threshold];
                    for (t, window) in cdf.windows(2).enumerate() {
                        assert!(
                            window[0] <= window[1],
                            "{label} decreasing={decreasing}: covariate {c} falls \
                             from {} to {} at threshold {t}",
                            window[0],
                            window[1]
                        );
                    }
                }

                let rows: Vec<(f64, f64)> = (0..n_covariate)
                    .map(|c| (fit.covariates[2 * c], fit.covariates[2 * c + 1]))
                    .collect();
                for i in 0..n_covariate {
                    for j in 0..n_covariate {
                        if i == j || !(rows[i].0 <= rows[j].0 && rows[i].1 <= rows[j].1) {
                            continue;
                        }
                        // For an increasing fit the CDF of the dominated covariate must
                        // dominate pointwise; a decreasing fit mirrors the roles.
                        let (upper, lower) = if decreasing { (j, i) } else { (i, j) };
                        for t in 0..n_threshold {
                            let high = fit.cdfs[upper * n_threshold + t];
                            let low = fit.cdfs[lower * n_threshold + t];
                            assert!(
                                high >= low,
                                "{label} decreasing={decreasing}: comparable pair \
                                 ({i}, {j}) inverted at threshold {t}: {high} < {low}"
                            );
                        }
                    }
                }
            }
        }
    }
}
