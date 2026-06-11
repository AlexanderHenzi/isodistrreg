use crate::functionals::{ClippingWrapper, KaplanMeier};
#[cfg(feature = "partial-order")]
use crate::partial_order;
#[cfg(feature = "partial-order")]
use crate::partial_order::CovariateGroups;
use crate::preprocessing::validate;
use crate::structures::{Decreasing, Increasing};
use crate::test::is_relative_eq_vec;
use crate::total_order::functionals::algorithm;
use crate::total_order::preprocessing::preprocess_censored;
use crate::total_order::preprocessing::preprocess_uncensored;
use crate::total_order::stochastic_dominance::{censored, uncensored};
use itertools::izip;

fn execute_test<const N: usize, const N_COVARIATE: usize, const N_THRESHOLD: usize>(
    x: [f64; N],
    y: [f64; N],
    weight: [f64; N],
    expected: [[f64; N_COVARIATE]; N_THRESHOLD],
) {
    let expected_flat: Vec<f64> = expected.iter().flatten().copied().collect();
    // f32-narrowed view for comparisons against the f32 algorithms (uncensored, censored
    // definition, censored fast). Done once, then reused.
    let expected_flat_f32: Vec<f32> = expected_flat.iter().map(|&v| v as f32).collect();

    {
        // Uncensored
        validate(x.chunks_exact(1), &y, None, Some(&weight)).unwrap();
        let context = preprocess_uncensored(&x, &y, &weight);
        let cdfs = uncensored::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.unique_responses.len(), N_THRESHOLD);
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );

        // Uncensored, reverse ordering
        let reversed_covariates = x.iter().map(|v| -v).collect::<Vec<_>>();
        validate(reversed_covariates.chunks_exact(1), &y, None, Some(&weight)).unwrap();
        let context_rev = preprocess_uncensored(&reversed_covariates, &y, &weight);
        let reversed_cdfs = uncensored::<Decreasing, _, _>(&context_rev, &crate::NoProgress);
        assert_eq!(context_rev.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context_rev.unique_responses.len(), N_THRESHOLD);
        let cdfs = reversed_cdfs
            .chunks_exact(N_COVARIATE)
            .map(|threshold| threshold.iter().rev())
            .flatten()
            .copied()
            .collect::<Vec<_>>();
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );
    }

    {
        // Censored
        let observed = [true; N];
        validate(x.chunks_exact(1), &y, Some(&observed), Some(&weight)).unwrap();

        // Definition
        let context = preprocess_censored(&x, &y, &observed, &weight).unwrap();
        let cdfs = censored::algorithm_definition::<Increasing, _, _>(&context);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.thresholds.len(), N_THRESHOLD);
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );

        // Definition, reversed
        let reversed_covariates = x.iter().map(|v| -v).collect::<Vec<_>>();
        let context_rev =
            preprocess_censored(&reversed_covariates, &y, &observed, &weight).unwrap();
        let reversed_cdfs = censored::algorithm_definition::<Decreasing, _, _>(&context_rev);
        assert_eq!(context_rev.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context_rev.thresholds.len(), N_THRESHOLD);
        let cdfs = reversed_cdfs
            .chunks_exact(N_COVARIATE)
            .map(|threshold| threshold.iter().rev().copied())
            .flatten()
            .collect::<Vec<_>>();
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );

        // Fast
        let context = preprocess_censored(&x, &y, &observed, &weight).unwrap();
        let cdfs = censored::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.thresholds.len(), N_THRESHOLD);
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );

        // Fast, reversed
        let context_rev =
            preprocess_censored(&reversed_covariates, &y, &observed, &weight).unwrap();
        let reversed_cdfs = censored::<Decreasing, _, _>(&context_rev, &crate::NoProgress);
        assert_eq!(context_rev.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context_rev.thresholds.len(), N_THRESHOLD);
        let cdfs = reversed_cdfs
            .chunks_exact(N_COVARIATE)
            .map(|threshold| threshold.iter().rev().copied())
            .flatten()
            .collect::<Vec<_>>();
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );
    }

    {
        // Censored, per threshold
        let observed = [true; N];
        let mut thresholds = Vec::from(y);
        thresholds.sort_unstable_by(|a, b| a.total_cmp(b));
        thresholds.dedup();
        assert_eq!(thresholds.len(), N_THRESHOLD);

        for (j, threshold) in thresholds.into_iter().enumerate() {
            let functional = ClippingWrapper::new(KaplanMeier::new(threshold));

            let result = algorithm::<Decreasing, _, _>(izip!(x, y, observed, weight), &functional);
            // `KaplanMeier` now produces `f32`, so the functional-PAVA result is `Vec<f32>` —
            // identical precision to what `fast::algorithm` exposes to users. Narrow the f64
            // reference array to f32 to compare like-with-like.
            let expected_f32: Vec<f32> = expected[j].iter().map(|&v| v as f32).collect();
            assert!(
                is_relative_eq_vec(&result, &expected_f32),
                "Result:   {:?}\nExpected: {:?}\n",
                result,
                expected[j],
            );
        }
    }

    #[cfg(feature = "partial-order")]
    {
        // Censored, partial order
        let observed = [true; N];
        let algorithm_context = partial_order::preprocess_censored(
            &x,
            &y,
            &observed,
            &weight,
            &CovariateGroups::empty(1),
        );
        assert_eq!(algorithm_context.n_threshold(), N_THRESHOLD);
        let output =
            partial_order::censored::<Increasing, _, _>(&algorithm_context, &crate::NoProgress);
        // The partial-order censored algorithm narrows its f64 accumulators to f32 at the return boundary.
        let expected_f32: Vec<f32> = expected_flat.iter().map(|&v| v as f32).collect();
        assert!(
            is_relative_eq_vec(&output.cdfs, &expected_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            output.cdfs,
            expected_flat,
        );
    }
}

#[test]
fn test_trivial_single_observation() {
    execute_test([1.0], [5.0], [2.0], [[1.0]]);
}

#[test]
fn test_trivial_single_covariate_2() {
    execute_test([2.0, 2.0], [3.0, 5.0], [1.0, 2.0], [[1.0 / 3.0], [1.0]]);
}

#[test]
fn test_trivial_single_covariate_3() {
    execute_test(
        [2.0, 2.0, 2.0],
        [3.0, 5.0, 6.0],
        [1.0, 2.0, 3.0],
        [[1.0 / 6.0], [1.0 / 2.0], [1.0]],
    );
}

#[test]
fn test_trivial_single_response_2() {
    execute_test([1.0, 2.0], [5.0, 5.0], [3.0, 4.0], [[1.0; 2]]);
}

#[test]
fn test_trivial_single_response_3() {
    execute_test(
        [1.0, 2.0, 3.0],
        [5.0, 5.0, 5.0],
        [3.0, 4.0, 5.0],
        [[1.0; 3]],
    );
}

#[test]
fn test_trivial_single_observation_duplicate() {
    execute_test([1.0, 1.0], [5.0, 5.0], [2.0, 8.1], [[1.0]]);
}

#[test]
fn test_monotone() {
    execute_test(
        [1.0, 2.0, 3.0],
        [1.0, 2.0, 3.0],
        [1.0, 2.0, 3.0],
        [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0; 3]],
    );
}

#[test]
fn test_not_monotone_3() {
    execute_test(
        [1.0, 2.0, 3.0],
        [4.0, 2.0, 3.0],
        [1.0; 3],
        [[0.5, 0.5, 0.0], [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0], [1.0; 3]],
    );
}

#[test]
fn test_not_monotone_4() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0],
        [2.0, 4.0, 3.0, 1.0],
        [1.0; 4],
        [
            [0.25, 0.25, 0.25, 0.25],
            [1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [1.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
            [1.0; 4],
        ],
    );
}

#[test]
fn test_duplicate_covariate() {
    execute_test(
        [1.0, 2.0, 2.0, 3.0],
        [1.0, 2.0, 3.0, 4.0],
        [1.0; 4],
        [[1.0, 0.0, 0.0], [1.0, 0.5, 0.0], [1.0, 1.0, 0.0], [1.0; 3]],
    );
}

#[test]
fn test_duplicate_response() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.5, 2.5, 4.0],
        [1.0; 4],
        [[1.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 0.0], [1.0; 4]],
    );
}

#[test]
fn test_duplicate_response_first() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 1.0, 2.0, 3.0],
        [1.0; 4],
        [
            [1.0, 1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
        ],
    );
}

#[test]
fn test_duplicate_response_last() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.0, 3.0, 3.0],
        [1.0; 4],
        [[1.0, 0.0, 0.0, 0.0], [1.0, 1.0, 0.0, 0.0], [1.0; 4]],
    );
}

#[test]
fn test_duplicate_response_triangle() {
    execute_test(
        [1.0, 3.0, 2.0],
        [3.0, 3.0, 5.0],
        [1.0; 3],
        [[1.0, 0.5, 0.5], [1.0; 3]],
    );
}

#[test]
fn test_weighted_duplicate() {
    execute_test(
        [0.0, 2.0, 1.0, 1.0, 3.0],
        [6.0, 7.0, 8.0, 8.0, 9.0],
        [1.0, 2.0, 1.0, 1.0, 1.0],
        [
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 0.5, 0.5, 0.0],
            [1.0, 1.0, 1.0, 0.0],
            [1.0; 4],
        ],
    );
}

#[test]
fn test_whether_restoring_partitions_can_be_accelerated() {
    execute_test(
        [0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0],
        [1.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0],
        [100.0, 9.0, 1.0, 8.0, 2.0, 7.0, 3.0, 6.0, 4.0],
        [
            [10.0 * (0.9 + 0.8 + 0.7 + 0.6) / 140.0; 5],
            [1.0, 0.9, 0.8, 0.7, 0.6],
            [1.0; 5],
        ],
    );
}

#[test]
fn test_5() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [3.0, 2.0, 5.0, 4.0, 1.0],
        [1.0; 5],
        [
            [0.2, 0.2, 0.2, 0.2, 0.2],
            [0.5, 0.5, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [1.0, 1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [1.0, 1.0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
            [1.0; 5],
        ],
    )
}

#[test]
fn test_zero_comparisons() {
    execute_test(
        [-1.0, -0.0, 0.0, 1.0],
        [-0.0, 5.0, 4.0, 0.0],
        [1.0; 4],
        [
            [1.0, 0.3333333333333333, 0.3333333333333333],
            [1.0, 0.6666666666666666, 0.6666666666666666],
            [1.0; 3],
        ],
    )
}
