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

/// Finalize the in-progress (last) observation: store its narrowed weight and fold it into
/// the current (last) covariate statistic.
///
/// A positive `W`-precision weight sum can narrow to exactly 0.0f32 (a `W = f64` weight
/// below half the smallest f32 subnormal). Such an observation carries no representable
/// mass in the f32 algorithm domain and is dropped here exactly like a zero-weight
/// observation (see the `fit()` weight contract), keeping the "every weight in the context
/// is positive" invariant intact.
fn commit_observation<Y, S>(
    observations: &mut Vec<Observation<usize, Y, S, f32>>,
    covariate_statistics: &mut [CovariateStatistic],
    finalized: f32,
) {
    if finalized > 0.0 {
        observations.last_mut().unwrap().weight = finalized;
        let statistic = covariate_statistics.last_mut().unwrap();
        statistic.weight += finalized;
        statistic.cumulative_weight += finalized;
    } else {
        observations.pop();
    }
}

/// Close the covariate whose observations have all been committed: if every one of its
/// observations was dropped (narrowed to 0.0f32), the covariate carries no mass and must not
/// enter the grid — exactly like a covariate whose observations all have zero weight.
fn close_covariate<X>(
    covariate_statistics: &mut Vec<CovariateStatistic>,
    unique_covariates: &mut Vec<X>,
) {
    if covariate_statistics.last().unwrap().weight == 0.0 {
        covariate_statistics.pop();
        unique_covariates.pop();
    }
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
///
/// Weights are narrowed to f32 here (see the `fit()` weight contract); positive weights
/// that narrow to 0.0f32 are dropped like zero weights ([`commit_observation`]), with the
/// covariate grid and threshold grid matching the post-drop subsample exactly.
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
            cumulative_weight: 0.0, // Accumulated as each observation's weight is committed
        });
        unique_covariates.push(x[index]);
    }

    for index in covariate_sorted_indices {
        let covariate_equal = x[index] == *unique_covariates.last().unwrap();
        let (response_equal, censoring_equal) = {
            let last_observation = observations.last().unwrap();
            (
                y[index] == last_observation.y,
                observed(index) == last_observation.observed,
            )
        };

        if covariate_equal && response_equal && censoring_equal {
            // At the same observation -> just accumulate observation weight in W precision
            last_w_accum = last_w_accum + weight[index];
        } else if covariate_equal {
            // New observation but same covariate -> finalize the previous observation's
            // weight (narrow to f32, dropping it if no mass survives the narrowing),
            // accumulate covariate weight from it, and start a new observation with a
            // fresh W accumulator.
            commit_observation(
                &mut observations,
                &mut covariate_statistics,
                last_w_accum.to_f32().unwrap(),
            );
            observations.push(Observation {
                x: unique_covariates.len() - 1, // we stay with the same covariate
                y: y[index],
                observed: observed(index),
                weight: 0.0, // placeholder; finalized from `last_w_accum`
            });
            last_w_accum = weight[index];
        } else {
            // New observation and new covariate -> finalize the previous observation's
            // weight, accumulate covariate weight from it, close the finished covariate
            // (dropping it if all its observations were dropped), push a new covariate
            // statistic chaining the cumulative weight, push a new observation, and seed
            // a fresh W accumulator.
            commit_observation(
                &mut observations,
                &mut covariate_statistics,
                last_w_accum.to_f32().unwrap(),
            );
            close_covariate(&mut covariate_statistics, &mut unique_covariates);
            let new_statistic = CovariateStatistic {
                weight: 0.0,
                cumulative_weight: covariate_statistics
                    .last()
                    .map_or(0.0, |statistic| statistic.cumulative_weight),
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
        }
    }
    // Finalize the last in-progress observation and close the last covariate.
    commit_observation(
        &mut observations,
        &mut covariate_statistics,
        last_w_accum.to_f32().unwrap(),
    );
    close_covariate(&mut covariate_statistics, &mut unique_covariates);
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
