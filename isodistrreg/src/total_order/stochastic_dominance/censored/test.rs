use crate::functionals::{ClippingWrapper, KaplanMeier};
#[cfg(feature = "partial-order")]
use crate::partial_order;
#[cfg(feature = "partial-order")]
use crate::partial_order::CovariateGroups;
use crate::preprocessing::validate;
use crate::structures::{Decreasing, Increasing};
use crate::test::is_relative_eq_vec;
use crate::total_order::functionals::algorithm;
use crate::total_order::stochastic_dominance::censored::preprocess;
use crate::total_order::stochastic_dominance::censored::{definition, fast};
use itertools::izip;

fn execute_test<const N: usize, const N_COVARIATE: usize, const N_THRESHOLD: usize>(
    x: [f64; N],
    y: [f64; N],
    weight: [f64; N],
    observed: [bool; N],
    expected: [[f64; N_COVARIATE]; N_THRESHOLD],
) {
    let expected_flat: Vec<f64> = expected.iter().flatten().copied().collect();
    // f32-narrowed view for comparisons against algorithms that compute in f32 (definition,
    // fast). Done once, then reused. The narrow happens through the literal value's f64
    // representation, so it matches what `as f32` produces from the algorithm's f64 inputs.
    let expected_flat_f32: Vec<f32> = expected_flat.iter().map(|&v| v as f32).collect();
    validate(x.chunks_exact(1), &y, Some(&observed), Some(&weight)).unwrap();

    fn unique_with_filter_sort<const N: usize>(vs: &[f64; N], filter: &[bool; N]) -> Vec<f64> {
        let mut filtered: Vec<_> = vs
            .iter()
            .zip(filter)
            .filter(|&(_, &b)| b)
            .map(|(&v, _)| v)
            .collect();
        filtered.sort_unstable_by(|l, r| l.total_cmp(r));
        filtered.dedup();
        filtered
    }

    // Definition
    {
        let context = preprocess(&x, &y, &observed, &weight).unwrap();
        let cdfs = definition::algorithm::<Increasing, _, _>(&context);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.thresholds.len(), N_THRESHOLD);
        assert_eq!(context.thresholds, unique_with_filter_sort(&y, &observed));
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat_f32,
        );
    }

    // Definition, reversed
    {
        let reversed_covariates = x.iter().map(|v| -v).collect::<Vec<_>>();
        let context = preprocess(&reversed_covariates, &y, &observed, &weight).unwrap();
        let reversed_cdfs = definition::algorithm::<Decreasing, _, _>(&context);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.thresholds.len(), N_THRESHOLD);
        assert_eq!(context.thresholds, unique_with_filter_sort(&y, &observed));
        if N_COVARIATE > 0 {
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
                expected_flat_f32,
            );
        }
    }

    // Fast
    {
        let context = preprocess(&x, &y, &observed, &weight).unwrap();
        let cdfs = fast::algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.thresholds, unique_with_filter_sort(&y, &observed));
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat_f32,
        );
    }

    // TODO: Fast, reversed
    // {
    //     let reversed_covariates = covariate.iter().map(|v| -v).collect::<Vec<_>>();
    //     let (reversed_cdfs, unique_covariates, thresholds) = fast::algorithm::<Decreasing, _, _>(&reversed_covariates, &response, &censoring, &weights);
    //     assert_eq!(unique_covariates.len(), N_COVARIATE);
    //     assert_eq!(thresholds.len(), N_THRESHOLD);
    //     let cdfs = reversed_cdfs
    //         .chunks_exact(N_COVARIATE)
    //         .map(|threshold| threshold.iter().rev())
    //         .flatten()
    //         .copied()
    //         .collect::<Vec<_>>();
    //     assert!(
    //         is_relative_eq_vec(&cdfs, &expected_flat),
    //         "Result:   {:?}\nExpected: {:?}\n", reversed_cdfs, expected_flat,
    //     );
    // }

    fn unique(vs: impl IntoIterator<Item = f64>) -> Vec<f64> {
        let mut sorted = vs.into_iter().collect::<Vec<_>>();
        sorted.sort_unstable_by(|l, r| l.total_cmp(r));
        sorted.dedup();
        sorted
    }

    // Per threshold, only execute if no covariates get dropped
    if N_COVARIATE == unique(x).len() {
        let thresholds = unique_with_filter_sort(&y, &observed);
        assert_eq!(thresholds.len(), N_THRESHOLD);

        if N_COVARIATE > 0 {
            for (threshold, expected) in thresholds.into_iter().zip(expected) {
                let functional = ClippingWrapper::new(KaplanMeier::new(threshold));

                let result =
                    algorithm::<Decreasing, _, _>(izip!(x, y, observed, weight,), &functional);
                // `KaplanMeier` now produces `f32`, so the functional-PAVA output is `Vec<f32>` —
                // same precision as the user-visible `fast` algorithm output. Narrow the f64
                // reference array to f32 to compare like-with-like.
                let expected_f32: Vec<f32> = expected.iter().map(|&v| v as f32).collect();
                assert!(
                    is_relative_eq_vec(&result, &expected_f32),
                    "Result:   {:?}\nExpected: {:?}\n",
                    result,
                    expected,
                );
            }
        }
    }

    #[cfg(feature = "partial-order")]
    {
        use std::iter::repeat_n;

        // Censored, partial order
        let dimensions = 3;
        let cov_repeated = x
            .iter()
            .map(|&c| repeat_n(c, dimensions))
            .flatten()
            .collect::<Vec<_>>();
        if !observed.iter().all(|&b| b) {
            let algorithm_context = partial_order::preprocess_censored(
                &cov_repeated,
                &y,
                &observed,
                &weight,
                &CovariateGroups::empty(dimensions),
            );
            let output =
                partial_order::censored::<Increasing, _, _>(&algorithm_context, &crate::NoProgress);
            // The partial-order censored algorithm narrows its f64 accumulators to f32 at the
            // return boundary; the test compares against f32-narrowed expectations.
            let expected_f32: Vec<f32> = expected_flat.iter().map(|&v| v as f32).collect();
            assert!(
                is_relative_eq_vec(&output.cdfs, &expected_f32),
                "Result:   {:?}\nExpected: {:?}\n",
                output.cdfs,
                expected_flat,
            );
        }
    }
}

#[test]
fn test_shortcuts_first_censored() {
    execute_test(
        [0.0, 1.0, 2.0, 3.0, 4.0],
        [0.0, 1.0, 2.0, 4.0, 3.0],
        [1.0; 5],
        [false, false, false, true, true],
        [[0.5, 0.5], [1.0, 1.0]],
    );
}

#[test]
fn test_monotone_censored() {
    execute_test(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        [1.0; 3],
        [true, false, true],
        [[1.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
    );
}

#[test]
fn test_monotone_censored_all_but_one() {
    execute_test(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        [1.0; 3],
        [false, true, false],
        [[1.0, 0.0]],
    );
}

#[test]
fn test_not_monotone_censored_3() {
    execute_test(
        [0.0, 1.0, 2.0],
        [4.0, 2.0, 3.0],
        [1.0; 3],
        [true, false, true],
        [[0.5, 0.5], [1.0, 1.0]],
    );
}

#[test]
fn test_monotone_partially_censored() {
    execute_test(
        [0.0, 1.0, 2.0, 3.0],
        [0.0, 1.0, 1.0, 2.0],
        [1.0; 4],
        [true, true, false, true],
        // There are two possible answers for this test at covariate 2.0, because the estimator
        // is not defined there. The alternative is:
        // [
        //     [1.0, 0.0, 0.0, 0.0],
        //     [1.0, 1.0, 1.0, 0.0],
        //     [1.0, 1.0, 1.0, 1.0],
        // ],
        [
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
        ],
    );
}

#[test]
fn test_monotone_partially_censored_reversed() {
    execute_test(
        [0.0, 1.0, 2.0, 3.0],
        [0.0, 1.0, 1.0, 2.0],
        [1.0; 4],
        [true, false, true, true],
        [
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 0.5, 0.5, 0.0],
            [1.0, 1.0, 1.0, 1.0],
        ],
    );
}

#[test]
fn test_monotone_partially_censored_duplicate_covariate() {
    execute_test(
        [0.0, 1.0, 1.0, 2.0],
        [0.0, 1.0, 1.0, 2.0],
        [1.0; 4],
        [true, true, false, true],
        [[1.0, 0.0, 0.0], [1.0, 0.5, 0.0], [1.0, 1.0, 1.0]],
    );
}

#[test]
fn test_monotone_partially_censored_duplicate_covariate_reversed() {
    execute_test(
        [0.0, 1.0, 1.0, 2.0],
        [0.0, 1.0, 1.0, 2.0],
        [1.0; 4],
        [true, false, true, true],
        [[1.0, 0.0, 0.0], [1.0, 0.5, 0.0], [1.0, 1.0, 1.0]],
    );
}

#[test]
fn test_not_monotone_first_censored() {
    execute_test(
        [184.0, 144.0],
        [0.2456797, 0.3201176],
        [1.0; 2],
        [true, false],
        [[0.5, 0.5]],
    );
}

#[test]
fn test_not_monotone_last_censored() {
    execute_test(
        [184.0, 144.0],
        [0.24, 0.32],
        [1.0; 2],
        [false, true],
        [[1.0]],
    );
}

#[test]
fn test_not_monotone_4() {
    execute_test(
        [6.0, 23.0, 44.0, 95.0],
        [0.4, 0.1, 0.6, 0.5],
        [1.0; 4],
        [false, true, true, true],
        [
            [0.5, 0.5, 0.0, 0.0],
            [0.5, 0.5, 0.5, 0.5],
            [1.0, 1.0, 1.0, 1.0],
        ],
    );
}

#[test]
fn test_4_1() {
    execute_test(
        [2.0, 1.0, 4.0, 3.0],
        [0.5, 4.0, 1.0, 0.0],
        [1.0; 4],
        [true, true, false, true],
        [
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
            [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 0.0],
            [1.0, 1.0, 1.0, 0.0],
        ],
    );
}

#[test]
///   1
///  0.5   0
///  0.5   0   1
/// 0.625  0  0.5  0
///   1
///  0.5  0
///  0.5  0   1
///  0.5  0  0.5  0
fn test_antitone_4() {
    execute_test(
        [2.0, 3.0, 1.0, 4.0],
        [3.0, 2.0, 4.0, 1.0],
        [1.0; 4],
        [true, false, false, true],
        [[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0.5, 0.5]],
    );
}

#[test]
fn test_7() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        [8.0, 4.0, 3.0, 2.0, 6.0, 5.0, 7.0],
        [1.0; 7],
        [true, true, true, true, false, false, true],
        [
            [0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0],
            [0.75, 0.75, 0.75, 0.75, 0.0, 0.0, 0.0],
            [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ],
    );
}

#[test]
fn test_9() {
    execute_test(
        [1.0, 3.0, 4.0, 6.0, 8.0, 5.0, 7.0, 9.0, 2.0],
        [7.0, 2.0, 6.0, 5.0, 9.0, 3.0, 8.0, 4.0, 1.0],
        [1.0; 9],
        [true, true, false, true, false, false, false, true, true],
        [
            [1.0 / 2.0, 1.0 / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [
                2.0 / 3.0,
                2.0 / 3.0,
                2.0 / 3.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            [
                2.0 / 3.0,
                2.0 / 3.0,
                2.0 / 3.0,
                0.2,
                0.2,
                0.2,
                0.2,
                0.2,
                0.2,
            ],
            [
                2.0 / 3.0,
                2.0 / 3.0,
                2.0 / 3.0,
                0.5,
                0.5,
                0.5,
                1.0 / 3.0,
                1.0 / 3.0,
                1.0 / 3.0,
            ],
            [
                1.0,
                1.0,
                1.0,
                0.5,
                0.5,
                0.5,
                1.0 / 3.0,
                1.0 / 3.0,
                1.0 / 3.0,
            ],
        ],
    );
}

#[test]
fn test_9_2() {
    execute_test(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [4.0, 9.0, 8.0, 5.0, 6.0, 7.0, 1.0, 3.0, 2.0],
        [1.0; 9],
        [false, false, false, true, false, true, true, false, true],
        [
            [
                1.0 / 7.0,
                1.0 / 7.0,
                1.0 / 7.0,
                1.0 / 7.0,
                1.0 / 7.0,
                1.0 / 7.0,
                1.0 / 7.0,
                0.0,
                0.0,
            ],
            [
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
                2.0 / 9.0,
            ],
            [
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
                17.0 / 45.0,
            ],
            [
                19.0 / 35.0,
                19.0 / 35.0,
                19.0 / 35.0,
                19.0 / 35.0,
                19.0 / 35.0,
                19.0 / 35.0,
                19.0 / 35.0,
                0.5,
                0.5,
            ],
        ],
    );
}

#[test]
fn test_18() {
    execute_test(
        [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
            17.0, 18.0,
        ],
        [
            16.0, 9.0, 17.0, 7.0, 10.0, 8.0, 2.0, 13.0, 15.0, 1.0, 14.0, 3.0, 11.0, 12.0, 4.0,
            18.0, 6.0, 5.0,
        ],
        [1.0; 18],
        [
            false, false, false, true, false, true, true, false, false, true, true, true, false,
            true, false, true, true, true,
        ],
        [
            [
                0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0,
            ],
            [
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0,
            ],
            [
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
            ],
            [
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.2, 0.2,
                0.2, 0.2, 0.2, 0.2,
            ],
            [2.0 / 7.0; 18],
            [29.0 / 84.0; 18],
            [
                3.0 / 7.0,
                3.0 / 7.0,
                3.0 / 7.0,
                3.0 / 7.0,
                3.0 / 7.0,
                3.0 / 7.0,
                3.0 / 7.0,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
                0.4,
            ],
            [24.0 / 49.0; 18],
            [67.0 / 112.0; 18],
            [1.0; 18],
        ],
    )
}

#[test]
fn test_fully_censored() {
    execute_test::<_, 0, _>(
        [
            -1.12910954,
            -0.82351673,
            -0.26303398,
            -0.06960184,
            1.99180191,
        ],
        [-2.2223270, -1.2314834, 0.3245255, 0.7515100, 1.0838672],
        [1.0; 5],
        [false; 5],
        [],
    )
}

#[test]
fn test_6x6() {
    execute_test(
        [
            0.3322026, 0.7690422, 0.5629895, -1.0083766, -0.3745809, -0.3724388,
        ],
        [
            -0.7637937, 0.2852616, 0.6007779, -1.2237571, -1.4237579, -1.6325940,
        ],
        [1.0; 6],
        [true, false, true, true, false, false],
        [
            [1.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.0],
        ],
    )
}
