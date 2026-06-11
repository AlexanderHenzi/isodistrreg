use crate::Float;
use crate::routines::argsort_indices_unstable_by;
use crate::structures::{Increasing, Observation};
use crate::total_order::structures::AlgorithmContext;
use crate::total_order::structures::CovariateStatistic;

pub fn preprocess_uncensored<X: Float, Y: Float, W: Float>(
    x: &[X],
    y: &[Y],
    weight: &[W],
) -> AlgorithmContext<X, Y, ()> {
    let mut context = preprocess(x, y, |_| (), weight);
    context.unique_responses = context
        .observations
        .first()
        .map(|o| o.y)
        .into_iter()
        .chain(
            context
                .observations
                .array_windows()
                .filter_map(|[left, right]| (right.y != left.y).then_some(right.y)),
        )
        .collect();
    context
}

/// Build the algorithm context: observations sorted and deduplicated, covariate grid and
/// per-covariate weight statistics.
///
/// Zero-weight observations are dropped here, before anything is aggregated: they carry no
/// statistical information (every estimator the kernels compute is invariant to adding
/// them), their response values create no thresholds, and a covariate whose observations
/// all have zero weight does not enter the grid. Downstream kernels divide by observation
/// and covariate weights and rely on every weight in the context being positive. With all
/// observations dropped, the returned context is empty (the empty fit).
pub fn preprocess<X: Float, Y: Float, S: Ord, W: Float, F: Fn(usize) -> S>(
    x: &[X],
    y: &[Y],
    observed: F,
    weight: &[W],
) -> AlgorithmContext<X, Y, S> {
    let n = x.len();

    // Determine the order by covariate and response, over the positive-weight
    // observations only — zero-weight ones are dropped before they are even compared.
    // `validate` has rejected negative and non-finite weights, so `> 0` is exactly
    // "nonzero" here.
    let order = argsort_indices_unstable_by::<Increasing, _>(
        |i, j| {
            x[i].total_cmp(&x[j])
                .then(y[i].total_cmp(&y[j]))
                .then(observed(i).cmp(&observed(j)).reverse())
        },
        (0..n).filter(|&i| weight[i] > W::zero()).collect(),
    );
    if order.is_empty() {
        return AlgorithmContext {
            observations: Vec::with_capacity(0),
            covariate_statistics: Vec::with_capacity(0),
            unique_responses: Vec::with_capacity(0),
            unique_covariates: Vec::with_capacity(0),
        };
    }

    // While copying the data in the sorted order and aggregating identical observations by weight,
    // we track the unique covariates we encounter.
    let mut observations: Vec<Observation<usize, Y, S, f32>> = Vec::with_capacity(n);
    let mut unique_covariates = Vec::with_capacity(n);
    let mut covariate_statistics: Vec<CovariateStatistic> = Vec::with_capacity(n);

    // Accumulator for the in-progress observation's weight in the caller's W precision.
    // Invariant: `last_w_accum` holds the running W-precision sum that will be committed
    // into `observations.last().weight` (narrowed to f32) when the in-progress observation
    // is finalized (i.e., when a new observation is pushed or at the end of the loop).
    // While the in-progress observation is current, its `weight: f32` field holds an
    // uninitialized placeholder (0.0) and must not be read.
    let mut last_w_accum: W;

    let mut covariate_sorted_indices = order.into_iter();
    {
        let index = covariate_sorted_indices.next().unwrap();
        observations.push(Observation {
            x: 0,
            y: y[index],
            observed: observed(index),
            weight: 0.0, // placeholder; finalized from `last_w_accum`
        });
        last_w_accum = weight[index];
        covariate_statistics.push(CovariateStatistic {
            weight: 0.0,
            cumulative_weight: 0.0, // We update the cumulative weight once duplicates are handled
        });
        unique_covariates.push(x[index]);
    }

    let mut current_statistic = &mut covariate_statistics[0];
    for index in covariate_sorted_indices {
        let covariate_equal = x[index] == *unique_covariates.last().unwrap();
        let last_observation = observations.last_mut().unwrap();
        let response_equal = y[index] == last_observation.y;
        let censoring_equal = observed(index) == last_observation.observed;

        if covariate_equal && response_equal && censoring_equal {
            // At the same observation -> just accumulate observation weight in W precision
            last_w_accum = last_w_accum + weight[index];
        } else if covariate_equal {
            // New observation but same covariate -> finalize the previous observation's
            // weight (narrow to f32), accumulate covariate weight from it, and start a new
            // observation with a fresh W accumulator.
            let finalized = last_w_accum.to_f32().unwrap();
            last_observation.weight = finalized;
            current_statistic.weight += finalized;
            observations.push(Observation {
                x: unique_covariates.len() - 1, // we stay with the same covariate
                y: y[index],
                observed: observed(index),
                weight: 0.0, // placeholder; finalized from `last_w_accum`
            });
            last_w_accum = weight[index];
        } else {
            // New observation and new covariate -> finalize the previous observation's
            // weight, accumulate covariate weight from it, collect the cumulative covariate
            // weight, push a new covariate statistic, push a new observation, and seed a
            // fresh W accumulator.
            let finalized = last_w_accum.to_f32().unwrap();
            last_observation.weight = finalized;
            current_statistic.weight += finalized;
            current_statistic.cumulative_weight += current_statistic.weight;
            let new_statistic = CovariateStatistic {
                weight: 0.0,
                cumulative_weight: current_statistic.cumulative_weight,
            };
            covariate_statistics.push(new_statistic);
            observations.push(Observation {
                x: unique_covariates.len(),
                y: y[index],
                observed: observed(index),
                weight: 0.0, // placeholder; finalized from `last_w_accum`
            });
            last_w_accum = weight[index];
            unique_covariates.push(x[index]);

            current_statistic = covariate_statistics.last_mut().unwrap();
        }
    }
    // Finalize the last in-progress observation.
    let finalized = last_w_accum.to_f32().unwrap();
    observations.last_mut().unwrap().weight = finalized;
    current_statistic.weight += finalized;
    current_statistic.cumulative_weight += current_statistic.weight;
    // Could shrink to fit `observations`, `unique_covariates` and `covariate_statistics` here if
    // duplicates are common, and we care about memory usage

    // Stable to keep each response threshold sorted by covariate - this is probably faster in PAVA
    observations.sort_by(|l, r| {
        l.y.total_cmp(&r.y)
            .then(l.observed.cmp(&r.observed).reverse())
    });

    AlgorithmContext {
        observations,
        covariate_statistics,
        unique_responses: Vec::with_capacity(0),
        unique_covariates,
    }
}

#[cfg(test)]
mod test {
    use crate::structures::Observation;
    use crate::total_order::preprocessing::preprocess_uncensored;
    use crate::total_order::structures::{AlgorithmContext, CovariateStatistic};

    #[test]
    fn test_trivial_single_observation() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[5.0], &[6.0], &[2.0]);
        assert_eq!(
            observations,
            vec![Observation {
                x: 0,
                y: 6.0,
                observed: (),
                weight: 2.0,
            },],
        );
        assert_eq!(
            covariate_statistics,
            vec![CovariateStatistic {
                weight: 2.0,
                cumulative_weight: 2.0,
            },],
        );
        assert_eq!(unique_responses, vec![6.0]);
        assert_eq!(unique_covariates, vec![5.0]);
    }

    #[test]
    fn test_trivial_single_covariate() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[5.0, 5.0], &[6.5, 6.0], &[1.0, 2.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 6.0,
                    observed: (),
                    weight: 2.0,
                },
                Observation {
                    x: 0,
                    y: 6.5,
                    observed: (),
                    weight: 1.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![CovariateStatistic {
                weight: 3.0,
                cumulative_weight: 3.0,
            }],
        );
        assert_eq!(unique_responses, vec![6.0, 6.5]);
        assert_eq!(unique_covariates, vec![5.0]);
    }

    #[test]
    fn test_trivial_single_response() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[5.0, 7.0], &[6.0, 6.0], &[1.0, 3.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 6.0,
                    observed: (),
                    weight: 1.0,
                },
                Observation {
                    x: 1,
                    y: 6.0,
                    observed: (),
                    weight: 3.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 1.0,
                },
                CovariateStatistic {
                    weight: 3.0,
                    cumulative_weight: 4.0,
                },
            ],
        );
        assert_eq!(unique_responses, vec![6.0]);
        assert_eq!(unique_covariates, vec![5.0, 7.0]);
    }

    #[test]
    fn test_monotone() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[2.0, 1.0, 3.0], &[2.0, 1.0, 3.0], &[1.0, 9.0, 1.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 1.0,
                    observed: (),
                    weight: 9.0,
                },
                Observation {
                    x: 1,
                    y: 2.0,
                    observed: (),
                    weight: 1.0,
                },
                Observation {
                    x: 2,
                    y: 3.0,
                    observed: (),
                    weight: 1.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![
                CovariateStatistic {
                    weight: 9.0,
                    cumulative_weight: 9.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 10.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 11.0,
                },
            ],
        );
        assert_eq!(unique_responses, vec![1.0, 2.0, 3.0]);
        assert_eq!(unique_covariates, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_not_monotone_3() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[2.0, 1.0, 3.0], &[3.0, 1.0, 2.0], &[1.0, 9.0, 1.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 1.0,
                    observed: (),
                    weight: 9.0,
                },
                Observation {
                    x: 2,
                    y: 2.0,
                    observed: (),
                    weight: 1.0,
                },
                Observation {
                    x: 1,
                    y: 3.0,
                    observed: (),
                    weight: 1.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![
                CovariateStatistic {
                    weight: 9.0,
                    cumulative_weight: 9.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 10.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 11.0,
                },
            ],
        );
        assert_eq!(unique_responses, vec![1.0, 2.0, 3.0]);
        assert_eq!(unique_covariates, vec![1.0, 2.0, 3.0]);
    }
}
