pub mod hazard_rate_order;

use crate::partial_order::algorithm::uncensored::hazard_rate_order::update_constraint_matrix;
use crate::partial_order::routines::derive_transitive_reduction;
use crate::partial_order::{
    AlgorithmContext, AlgorithmOutput, Config, OrderingInfo, QualityIndicators,
};
use crate::progress::ProgressTracker;
use crate::routines::transpose;
use crate::structures::{Direction, Increasing};
use crate::total_order::tonic_regression_pre_sorted;
use std::borrow::Cow;
use std::iter::repeat_n;

#[must_use]
pub fn algorithm<D: Direction, const HRO: bool>(
    context: &AlgorithmContext,
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
                    Some(*acc)
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

    // Build A `(constraints.len() x n_covariate)` matrix in CSC. We do an antitonic
    // (non-increasing) regression, so we reverse the order.
    let mut constraint_matrix =
        build_order_constraints::<D::REVERSE, HRO>(&constraint_edges, context.n_covariate());

    // Build diagonal P = diag(w_groups) (upper triangular only)
    let weight_matrix = build_diagonal_matrix(&context.x_weight);

    // State variables maintained as we iterate over the data
    let mut data_index = 0;
    // Construct the linear cost vector `q`. It is the negative of the element-wise product of
    // (grouped) weights and mean response. Together with the diagonal P, it represents each
    // component of the loss like (with p_i share of weight below current threshold):
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
    // Under HRO, we work in the reverse direction because we model the survival function.
    let mut q = if HRO {
        context.x_weight.iter().map(|w| -w).collect()
    } else {
        vec![0.0; context.n_covariate()]
    };

    // We collect these results each iteration
    let mut cdfs = Vec::with_capacity(context.n_threshold() * context.n_covariate());
    let mut iter_limit_hit_count = 0;

    // First iteration
    let mut current_threshold = None;
    while data_index < context.n() {
        match current_threshold {
            None => current_threshold = Some(context.y[data_index]),
            Some(existing_value) => {
                if context.y[data_index] > existing_value {
                    break;
                }
            }
        }
        let covariate_index = context.x_indices[data_index];
        if HRO {
            q[covariate_index] += context.weight[data_index];
        } else {
            q[covariate_index] -= context.weight[data_index];
        }
        data_index += 1;
    }

    let mut problem = osqp::Problem::new(
        weight_matrix,
        &q,
        // TODO: This cloning only happens for the HRO case, is it optimized away for the SD case?
        constraint_matrix.clone(),
        &vec![0.0; constraint_edges.len()],
        // TODO: The upper bound 1.0 is redundant, but in v1.0.0 osqp introduced a termination check
        //  via a duality-gap test which can't be computed without an upper bound. Once this
        //  termination criterion can be turned off via the rust api, this upper bound can be set to
        //  f64::INFINITY also (although it should be tested whether this actually improves the
        //  performance). The osqp solution then (still) should be clamped up to 1.0 to counteract
        //  numerical errors as a post-processing step regardless.
        &vec![1.0; constraint_edges.len()],
        &config.osqp_settings,
    )
    .expect("Failed to setup OSQP problem");

    // Initial solve
    problem.warm_start_x(&vec![if HRO { 1.0 } else { 0.0 }; context.n_covariate()]);
    let status = problem.solve();
    if let osqp::Status::MaxIterationsReached(_) = status {
        iter_limit_hit_count += 1;
    }
    let solution = status.solution().expect("Need OSQP to find a solution");
    let old_length = cdfs.len();
    cdfs.extend(solution.x().iter().map(|v| v.clamp(0.0, 1.0)));
    progress.increment();
    let mut primal_variable = &cdfs[old_length..];

    let last_threshold = *context.thresholds.last().unwrap();
    loop {
        // Bounds check not needed, because last threshold is not yet done
        let active_threshold = context.y[data_index];
        if active_threshold == last_threshold {
            break;
        }

        while context.y[data_index] == active_threshold {
            let covariate_index = context.x_indices[data_index];
            if HRO {
                q[covariate_index] += context.weight[data_index];
            } else {
                q[covariate_index] -= context.weight[data_index];
            }
            data_index += 1;
        }

        // Warm-started remaining non-trivial solves. OSQP warm-starts by default, but we prefer to
        // restart at the value clamped between 0 and 1.
        problem.warm_start_x(primal_variable);
        // Pass the updated linear cost
        problem.update_lin_cost(&q);
        if HRO {
            update_constraint_matrix::<D>(&mut constraint_matrix, primal_variable);
            problem.update_A(constraint_matrix.clone());
        }
        let status = problem.solve();
        if let osqp::Status::MaxIterationsReached(_) = status {
            iter_limit_hit_count += 1;
        }
        let solution = status.solution().expect("Need OSQP to find a solution");
        let old_length = cdfs.len();
        cdfs.extend(solution.x().iter().map(|v| v.clamp(0.0, 1.0)));
        progress.increment();
        primal_variable = &cdfs[old_length..];
    }

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
        let increasing = tonic_regression_pre_sorted::<Increasing>(cdf.iter().map(|&v| (v, 1.0)));
        for (initial, cleaned) in cdf.iter_mut().zip(increasing) {
            *initial = cleaned;
        }
    }

    let convergence_fraction = 1.0 - iter_limit_hit_count as f64 / context.n_threshold() as f64;

    AlgorithmOutput {
        cdfs,
        ordering_info: OrderingInfo::from_edges(constraint_edges, context.n_covariate()),
        quality_indicators: QualityIndicators {
            precision,
            convergence_fraction,
        },
    }
}

/// Build A as CSC from constraints (x_i <= x_j).
///
/// We create one row per edge with -1 at i and +1 at j for the constraint a^T x = x_j - x_i >= 0.
/// If HRO, the signs of these coefficients get flipped, because then we're fitting a survival curve
/// and not a CDF.
///
/// A has shape m x n, with m = number of edges and n = number of unique covariate rows.
pub fn build_order_constraints<D: Direction, const HRO: bool>(
    constraints: &[(usize, usize)],
    n_covariate: usize,
) -> osqp::CscMatrix<'_> {
    let m = constraints.len();
    let n = n_covariate;

    // Count nonzeros per column (each constraint contributes to two columns).
    let mut col_counts = vec![0usize; n];
    for &(i, j) in constraints {
        col_counts[i] += 1;
        col_counts[j] += 1;
    }

    // Build indptr via exclusive prefix sum of counts.
    let mut indptr = Vec::with_capacity(n + 1);
    indptr.push(0);
    for c in 0..n {
        let next = indptr[c] + col_counts[c];
        indptr.push(next);
    }
    let nnz = *indptr.last().unwrap_or(&0);

    // Allocate CSC storage.
    let mut indices = vec![0usize; nnz];
    let mut data = vec![0.0f64; nnz];

    // Running insertion pointers per column start at indptr[c].
    let mut next = indptr[..n].to_vec();

    // Fill columns in one pass. Row indices per column are increasing by construction.
    for (r, &(i, j)) in constraints.iter().enumerate() {
        let p_i = next[i];
        indices[p_i] = r;
        data[p_i] = if D::IS_INCREASING != HRO { -1.0 } else { 1.0 };
        next[i] += 1;

        let p_j = next[j];
        indices[p_j] = r;
        data[p_j] = if D::IS_INCREASING != HRO { 1.0 } else { -1.0 };
        next[j] += 1;
    }

    osqp::CscMatrix {
        nrows: m,
        ncols: n,
        indptr: Cow::Owned(indptr),
        indices: Cow::Owned(indices),
        data: Cow::Owned(data),
    }
}

/// Build a diagonal `osqp::CscMatrix` (no such method is available in the OSQP interface,
/// surprisingly).
fn build_diagonal_matrix(diag: &[f64]) -> osqp::CscMatrix<'_> {
    let n = diag.len();
    osqp::CscMatrix {
        nrows: n,
        ncols: n,
        indptr: Cow::Owned((0..=n).collect()),
        indices: Cow::Owned((0..n).collect()),
        data: Cow::Owned(diag.to_vec()),
    }
}

#[cfg(test)]
mod test {
    use crate::IsotonicDistributionalRegressionFit;
    use crate::partial_order::structures::Fit;
    use crate::partial_order::{Config, CovariateGroups, Csr, OrderingInfo};
    use crate::structures::StochasticOrder;

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
            osqp_settings: osqp::Settings::default()
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
            osqp_settings: osqp::Settings::default()
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
                osqp_settings: osqp::Settings::default()
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
