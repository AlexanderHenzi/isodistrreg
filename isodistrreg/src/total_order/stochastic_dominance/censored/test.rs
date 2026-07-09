use crate::functionals::{ClippingWrapper, KaplanMeier};
#[cfg(feature = "partial-order")]
use crate::partial_order;
#[cfg(feature = "partial-order")]
use crate::partial_order::CovariateGroups;
use crate::preprocessing::validate;
use crate::structures::{Decreasing, Increasing};
use crate::test::is_relative_eq_vec;
use crate::total_order::functionals::algorithm;
use crate::total_order::preprocessing::preprocess_censored as preprocess;
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

/// An event directly following a censored observation at the same covariate keeps
/// its own threshold and event mass. Thresholds [0, 2]; the x=2 cell dies at t=2.
#[test]
fn test_event_after_censored_same_covariate() {
    execute_test(
        [1.0, 2.0, 2.0],
        [0.0, 1.0, 2.0],
        [1.0; 3],
        [true, false, true],
        [[1.0, 0.0], [1.0, 1.0]],
    );
}

/// Kaplan-Meier factors are weight ratios, so the fit is invariant under scaling
/// all weights by a positive constant — including scales far below 1.
#[test]
fn test_weight_scale_invariance() {
    let expected = [[0.5, 0.0], [1.0, 0.0], [1.0, 1.0]];
    for scale in [1.0, 1e-6] {
        execute_test(
            [1.0, 1.0, 2.0, 2.0],
            [1.0, 2.0, 1.5, 3.0],
            [scale; 4],
            [true, true, false, true],
            expected,
        );
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

/// Regression test for double-applied Kaplan-Meier factors on tied uncensored responses.
///
/// The two uncensored observations tied at y=3 (covariates 2 and 4) form one run. While the
/// first of them is processed, `update_value`'s catch-up walk reads ahead and absorbs the
/// *entire* run into every cell it touches; the bug was recording `cold.count = data_index + 1`
/// (an index *inside* the run) instead of the run's end, so cells touched again while
/// processing the second tied observation re-applied its tail factors. On this instance the
/// broken bookkeeping yielded `[1.0, 0.5, 0.5]` instead of `[1.0, 1.0, 1.0]` at the second
/// threshold.
#[test]
fn test_tied_uncensored_responses() {
    execute_test(
        [4.0, 2.0, 4.0, 3.0],
        [1.0, 3.0, 3.0, 1.0],
        [1.0; 4],
        [true, true, true, false],
        [[0.25, 0.25, 0.25], [1.0, 1.0, 1.0]],
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

/// Differential tests: `fast` must agree with `definition` (the executable specification)
/// on generated instances. Comparison uses a small absolute tolerance — the two
/// implementations evaluate the same quantities with different operation orderings (and
/// `definition` accumulates in f64), so agreement is only expected up to f32 rounding. The
/// observed maximum deviation over hundreds of thousands of generated instances is ~1.5e-6
/// (about a dozen f32 ulps, on the largest continuous-data cases); the historical bugs these
/// tests guard against produced errors of 1e-2 and larger.
mod differential {
    use crate::structures::Increasing;
    use crate::total_order::preprocessing::preprocess_censored as preprocess;
    use crate::total_order::stochastic_dominance::censored::{definition, fast};
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    const TOLERANCE: f32 = 1e-5;

    /// Observation weights away from both 0 and uniformity.
    fn random_weight(rng: &mut StdRng) -> f64 {
        rng.random_range(0.5..2.0)
    }

    struct Instance {
        x: Vec<f64>,
        y: Vec<f64>,
        observed: Vec<bool>,
        weights: Vec<f64>,
    }

    impl Instance {
        /// Covariates and responses drawn from small integer grids: tiny level counts
        /// maximize tied responses, duplicate covariates, and censoring interactions per
        /// instance.
        fn integer_grid(
            rng: &mut StdRng,
            n: usize,
            covariate_levels: u32,
            response_levels: u32,
            censored_percent: u32,
            random_weights: bool,
        ) -> Self {
            Self {
                x: (0..n)
                    .map(|_| rng.random_range(1..=covariate_levels) as f64)
                    .collect(),
                y: (0..n)
                    .map(|_| rng.random_range(1..=response_levels) as f64)
                    .collect(),
                observed: (0..n)
                    .map(|_| rng.random_range(0..100) >= censored_percent)
                    .collect(),
                weights: (0..n)
                    .map(|_| {
                        if random_weights {
                            random_weight(rng)
                        } else {
                            1.0
                        }
                    })
                    .collect(),
            }
        }

        /// Survival-shaped data: the response is a monotone signal in the covariate plus
        /// independent noise (`rho` = signal share), right-censoring times are drawn below
        /// the event time, and covariate/response are optionally quantized onto discrete
        /// grids.
        fn survival(
            rng: &mut StdRng,
            n: usize,
            covariate_levels: Option<u32>,
            response_levels: Option<u32>,
            rho: f64,
            censoring: f64,
            random_weights: bool,
        ) -> Self {
            let quantize = |value: f64, levels: Option<u32>| match levels {
                Some(levels) => (value * levels as f64).floor() / levels as f64,
                None => value,
            };
            let mut instance = Instance {
                x: Vec::with_capacity(n),
                y: Vec::with_capacity(n),
                observed: Vec::with_capacity(n),
                weights: Vec::with_capacity(n),
            };
            for _ in 0..n {
                let x = quantize(rng.random(), covariate_levels);
                let event_time = rho * x + (1.0 - rho) * rng.random::<f64>();
                let censored = rng.random_bool(censoring);
                let y = if censored {
                    event_time * rng.random::<f64>()
                } else {
                    event_time
                };
                instance.x.push(x);
                instance.y.push(quantize(y, response_levels));
                instance.observed.push(!censored);
                instance.weights.push(if random_weights {
                    random_weight(rng)
                } else {
                    1.0
                });
            }
            instance
        }

        /// Run all implementations and compare; `Err` describes the worst difference and
        /// the full instance.
        fn check(&self) -> Result<(), String> {
            self.check_full(true)
        }

        /// Compare `definition` against `fast`; with `run_partial` (and the
        /// `partial-order` feature) additionally route the instance through the
        /// partial-order censored solver — the 1-D covariates replicated to three
        /// identical dimensions form a chain in the componentwise order — and through
        /// the literal min-max/RKM specification
        /// (`partial_order::algorithm::definition`), all against the same expected
        /// values.
        #[cfg_attr(not(feature = "partial-order"), allow(unused_variables))]
        fn check_full(&self, run_partial: bool) -> Result<(), String> {
            let context = preprocess(&self.x, &self.y, &self.observed, &self.weights).unwrap();
            let expected = definition::algorithm::<Increasing, _, _>(&context);
            let actual = fast::algorithm::<Increasing, _, _>(&context, &crate::NoProgress);

            assert_eq!(expected.len(), actual.len());
            let worst = expected
                .iter()
                .zip(&actual)
                .map(|(e, a)| (e - a).abs())
                .fold(0.0f32, f32::max);
            if worst > TOLERANCE {
                return Err(format!(
                    "max |definition - fast| = {worst:e}\n  \
                     x = {:?}\n  y = {:?}\n  observed = {:?}\n  weights = {:?}\n  \
                     definition = {expected:?}\n  fast = {actual:?}",
                    self.x, self.y, self.observed, self.weights,
                ));
            }

            #[cfg(feature = "partial-order")]
            if run_partial {
                self.check_partial(&expected)?;
            }
            Ok(())
        }

        /// Chain correspondence of the general partial-order algorithms with the
        /// total-order `definition`.
        #[cfg(feature = "partial-order")]
        fn check_partial(&self, expected: &[f32]) -> Result<(), String> {
            use crate::partial_order::{self, CovariateGroups};
            use std::iter::repeat_n;

            let dimensions = 3;
            let cov_repeated: Vec<f64> = self
                .x
                .iter()
                .flat_map(|&c| repeat_n(c, dimensions))
                .collect();
            let context = partial_order::preprocess_censored(
                &cov_repeated,
                &self.y,
                &self.observed,
                &self.weights,
                &CovariateGroups::empty(dimensions),
            );

            let solver =
                partial_order::censored::<Increasing, _, _>(&context, &crate::NoProgress).cdfs;
            self.compare_partial("partial-order solver", &solver, expected)?;

            // The literal specification enumerates (lower set, upper set) pairs:
            // polynomial on chains, but the per-subset partition enumeration is
            // exponential in the interval length, hence the covariate-count cap.
            if context.n_covariate() <= 12 {
                let spec = partial_order::censored_definition::<Increasing, _, _>(&context);
                self.compare_partial("min-max specification", &spec, expected)?;
            }
            Ok(())
        }

        #[cfg(feature = "partial-order")]
        fn compare_partial(
            &self,
            label: &str,
            actual: &[f32],
            expected: &[f32],
        ) -> Result<(), String> {
            let instance = format!(
                "x = {:?}\n  y = {:?}\n  observed = {:?}\n  weights = {:?}",
                self.x, self.y, self.observed, self.weights,
            );
            if actual.len() != expected.len() {
                return Err(format!(
                    "{label}: output length {} != total-order definition length {}\n  {instance}",
                    actual.len(),
                    expected.len(),
                ));
            }
            let worst = expected
                .iter()
                .zip(actual)
                .map(|(e, a)| (e - a).abs())
                .fold(0.0f32, f32::max);
            if worst <= TOLERANCE {
                Ok(())
            } else {
                Err(format!(
                    "max |definition - {label}| = {worst:e}\n  {instance}\n  \
                     definition = {expected:?}\n  {label} = {actual:?}",
                ))
            }
        }
    }

    fn assert_all_ok(failures: Vec<String>) {
        if let Some(first) = failures.first() {
            panic!("{} failing instances, first:\n{first}", failures.len());
        }
    }

    /// Tie-heavy instances large enough that every threshold carries a long run of tied
    /// uncensored responses, exercising the batched tie-run update against `definition`
    /// (the per-arrival and batched schedules must evaluate the same per-threshold
    /// estimator).
    #[test]
    fn tied_response_runs() {
        let mut failures = Vec::new();
        let mut seed = 4000u64;
        for n in [120, 250] {
            for response_levels in [3, 12, 40] {
                for censoring in [0.0, 0.5, 0.8] {
                    for covariate_levels in [None, Some(25)] {
                        for rho in [0.0, 0.5] {
                            seed += 1;
                            let instance = Instance::survival(
                                &mut StdRng::seed_from_u64(seed),
                                n,
                                covariate_levels,
                                Some(response_levels),
                                rho,
                                censoring,
                                seed.is_multiple_of(2),
                            );
                            failures.extend(
                                instance
                                    .check_full(false)
                                    .err()
                                    .map(|e| format!("seed {seed}: {e}")),
                            );
                        }
                    }
                }
            }
        }
        assert_all_ok(failures);
    }

    /// Factor grid matching the benchmark axes: continuous/discrete covariate and response,
    /// strong/weak/no X-Y dependence, low/high censoring share, several sizes.
    #[test]
    fn factor_grid() {
        let mut failures = Vec::new();
        let mut seed = 1u64;
        for covariate_levels in [None, Some(7)] {
            for response_levels in [None, Some(9)] {
                for rho in [0.9, 0.3, 0.0] {
                    for censoring in [0.3, 0.7] {
                        for n in [1, 2, 3, 5, 17, 60] {
                            seed += 1;
                            let instance = Instance::survival(
                                &mut StdRng::seed_from_u64(seed),
                                n,
                                covariate_levels,
                                response_levels,
                                rho,
                                censoring,
                                seed.is_multiple_of(2),
                            );
                            failures.extend(
                                instance.check().err().map(|e| format!("seed {seed}: {e}")),
                            );
                        }
                    }
                }
            }
        }
        assert_all_ok(failures);
    }

    /// Dense micro-instance sweep. This is the search that found
    /// `test_tied_uncensored_responses` and `clip_order_instance` (the factor grid ran at
    /// sizes where those patterns slipped through).
    #[test]
    fn micro_instances() {
        let mut rng = StdRng::seed_from_u64(0x5eed);
        let mut failures = Vec::new();
        for n in 3..=8 {
            for attempt in 0usize..2000 {
                let instance =
                    Instance::integer_grid(&mut rng, n, 4, 3, 40, !attempt.is_multiple_of(2));
                failures.extend(instance.check().err());
            }
        }
        assert_all_ok(failures);
    }

    /// Deep brute-force sweep: ~210 000 deterministic instances over n up to 40, 2-10
    /// covariate levels, 2-8 response levels, 20-70% censoring, equal and random weights.
    /// The per-instance cost is roughly linear in n, so attempts scale as 1/n to spend
    /// comparable time on every size; the whole sweep stays within seconds even in debug
    /// builds.
    #[test]
    fn deep_search() {
        let mut rng = StdRng::seed_from_u64(0xfeedbeef);
        let mut failures = Vec::new();
        for n in [3, 4, 5, 6, 7, 8, 10, 12, 18, 25, 40] {
            for attempt in 0..140_000 / n {
                let instance = Instance::integer_grid(
                    &mut rng,
                    n,
                    [2, 4, 6, 10][attempt % 4],
                    [2, 3, 5, 8][(attempt / 4) % 4],
                    [20, 40, 70][(attempt / 16) % 3],
                    attempt.is_multiple_of(2),
                );
                // The partial-order legs run on a 1-in-13 subsample to keep this
                // sweep within its in-debug seconds budget; 13 is coprime to the
                // %4 / %16 / %3 parameter cycling (period 48), so the subsample
                // still covers every parameter combination. Full partial-order
                // coverage at smaller scale comes from `factor_grid` and
                // `micro_instances`.
                failures.extend(instance.check_full(attempt % 13 == 0).err());
            }
        }
        assert_all_ok(failures);
    }

    /// Regression instance for the clip-order divergence. `definition` originally clipped
    /// in (r asc, s asc, k asc) order, where the column-side inputs `(k+1, s)` were not yet
    /// clipped at the current threshold, while `fast` clips against fully-clipped values on
    /// both sides; on this instance the two orders gave visibly different results (CDF
    /// 0.6549 vs 0.6308 at the second threshold). Resolved by defining clipping in
    /// span-ascending order, so clips only ever use already-clipped values — matching
    /// `fast`.
    #[test]
    fn clip_order_instance() {
        let instance = Instance {
            x: vec![1.0, 2.0, 2.0, 1.0, 3.0, 1.0, 4.0, 3.0],
            y: vec![1.0, 1.0, 3.0, 1.0, 2.0, 1.0, 1.0, 1.0],
            observed: vec![false, true, false, false, true, true, false, false],
            weights: vec![
                0.5948106348514557,
                1.312944918870926,
                0.6167095303535461,
                0.825943648815155,
                0.5218165516853333,
                1.0982480347156525,
                0.674931526184082,
                1.6719395518302917,
            ],
        };
        instance.check().unwrap();
    }
}

/// Regression test for exact-1.0 pinning at *intermediate* thresholds
/// (`CompletionIndex` interleaved into the fast algorithm; the deleted
/// `pin_completed_mass` post-pass only repaired the final threshold row).
///
/// Instance (sorted covariates 1..5): cov 1: event y=2 (w=1.7); cov 2: event y=1
/// (w=0.3); cov 3: censored y=4; cov 4: event y=5; cov 5: event y=3. The sub-CDFs at
/// covariates 1 and 2 complete at threshold y=2 — every interval Kaplan-Meier cell over
/// `{1}`, `{2}`, or `[1, 2]` ends in an event there, so the CDF is mathematically
/// exactly 1.0 from the second threshold on. The drift-prone weights make the f32
/// running sums land off by ulps: without the pin, the row at the *intermediate*
/// threshold y=3 came out as `[1.0, 0.9999998, ...]` (the bridged diagonal cell's
/// BIT-reconstructed consumed weight differs from the covariate weight by an ulp), so
/// the column at covariate 2 DECREASED along the thresholds.
#[test]
fn pins_exact_one_at_intermediate_threshold() {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [2.0, 1.0, 4.0, 5.0, 3.0];
    let observed = [true, true, false, true, true];
    let w = [1.7, 0.3, 1.0, 1.0, 1.0];

    let context = preprocess(&x, &y, &observed, &w).unwrap();
    let n_covariate = context.n_covariate();
    let n_threshold = context.n_threshold();
    assert_eq!((n_covariate, n_threshold), (5, 4));
    assert_eq!(context.thresholds, vec![1.0, 2.0, 3.0, 5.0]);

    let cdfs = fast::algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
    for covariate in 0..2 {
        // Exactly 1.0 (bitwise) from the second threshold on — including the
        // intermediate rows.
        for threshold in 1..n_threshold {
            assert_eq!(
                cdfs[threshold * n_covariate + covariate],
                1.0,
                "covariate {covariate}, threshold index {threshold}: completed sub-CDF \
                 must be exactly 1.0; cdfs = {cdfs:?}",
            );
        }
    }
    // Every column must be nondecreasing along the thresholds.
    for covariate in 0..n_covariate {
        let column: Vec<f32> = (0..n_threshold)
            .map(|threshold| cdfs[threshold * n_covariate + covariate])
            .collect();
        assert!(
            column.is_sorted(),
            "column at covariate index {covariate} decreases in the threshold: {column:?}",
        );
    }
}

/// The companion invariant across covariates: for an increasing fit every threshold
/// row is nonincreasing along the sorted covariates. Cells whose sub-CDF completes
/// must be pinned to exactly 1.0 even when their *neighbors* arrive at 1.0 through a
/// different route (e.g. the clamp of a slightly-negative survival) — otherwise a
/// drifted 0.99999994 between two exact 1.0s breaks the stochastic ordering by an ulp.
/// Here the first three covariates complete at the first threshold.
#[test]
fn pins_exact_one_preserves_covariate_ordering() {
    let x = [
        0.04097593542105893,
        0.05876877680199277,
        0.019171183910507095,
        0.16271247105566322,
        0.7878839749351719,
    ];
    let y = [
        0.12012876644649723,
        0.10460697733118048,
        0.10047503132764332,
        0.19127475492081214,
        0.7837285927709854,
    ];
    let observed = [true, true, true, false, true];
    let w = [1.4884429840234832, 1.0, 1.0, 1.0, 1.0];

    let context = preprocess(&x, &y, &observed, &w).unwrap();
    let n_covariate = context.n_covariate();
    let cdfs = fast::algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
    for (threshold, row) in cdfs.chunks_exact(n_covariate).enumerate() {
        assert!(
            row.iter().rev().is_sorted(),
            "threshold index {threshold}: row must be nonincreasing along the \
             covariates, got {row:?}",
        );
    }
    // At the third threshold (y = 0.12012...) all three left sub-CDFs have completed.
    assert_eq!(
        &cdfs[2 * n_covariate..3 * n_covariate],
        &[1.0, 1.0, 1.0, 0.0, 0.0]
    );
}

/// `preprocess` merges an *uncensored* observation into an immediately preceding
/// *censored* observation at the same covariate (`is_following_censored` in
/// preprocessing.rs), converting the event into censoring mass at a lower response and
/// silently dropping its response value from `thresholds`.
///
/// Derivation (single covariate, so S-IDR must reduce to the standard weighted
/// Kaplan-Meier estimator; there is no isotonicity constraint to satisfy):
///   data at the one covariate: event@0 (w=1), censored@1 (w=1), event@2 (w=1)
///   thresholds = unique uncensored responses = [0, 2]
///   (also per the `preprocess` doc: "Thresholds contain all unique response values,
///    discarding thresholds containing only censored observations" — y=2 is uncensored)
///   K-M: S(0) = 1 - 1/3 = 2/3            => CDF(0) = 1/3
///        S(2) = 2/3 * (1 - 1/1) = 0      => CDF(2) = 1
///   (risk set at time 2 is only the event itself: the censored point left at t=1)
///   expected output: thresholds [0, 2], cdfs [1/3, 1]
#[test]
fn preprocess_merges_event_into_censored_same_covariate() {
    let context = preprocess(
        &[1.0, 1.0, 1.0],
        &[0.0, 1.0, 2.0],
        &[true, false, true],
        &[1.0; 3],
    )
    .unwrap();
    let cdfs = fast::algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
    assert_eq!(
        context.thresholds,
        vec![0.0, 2.0],
        "uncensored response 2.0 must be a threshold"
    );
    assert_eq!(cdfs.len(), 2);
    assert!(
        (cdfs[0] - 1.0 / 3.0).abs() < 1e-6,
        "CDF(0) = 1/3, got {cdfs:?}"
    );
    assert!((cdfs[1] - 1.0).abs() < 1e-6, "CDF(2) = 1, got {cdfs:?}");
}
