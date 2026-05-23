use crate::Float;
use crate::error::Error;
use crate::preprocessing::validate;
use crate::routines::argsort_unstable_by;
use crate::structures::{Increasing, Observation};
use crate::total_order::stochastic_dominance::censored::structures::CensoredSdContext;
use crate::total_order::structures::CovariateStatistic;
use itertools::Itertools;

/// Preprocess data such that the algorithm can be run.
///
/// Thresholds contain all unique response values, discarding thresholds containing only censored
/// observations.
pub fn preprocess<X: Float, Y: Float, W: Float>(
    x: &[X],
    y: &[Y],
    observed: &[bool],
    weights: &[W],
) -> Result<CensoredSdContext<X, Y>, Error> {
    let n = validate(x.chunks_exact(1), y, Some(observed), Some(weights))?;

    let (observations_response_sorted, thresholds) = {
        let response_order = argsort_unstable_by::<Increasing, _>(
            |i, j| {
                y[i].total_cmp(&y[j])
                    .then(observed[i].cmp(&observed[j]).reverse())
                    .then(x[i].total_cmp(&x[j]))
            },
            n,
        );

        // Discard censored observations not greater than or equal to any uncensored observation
        let first_uncensored = response_order
            .iter()
            .find_position(|&&i| observed[i])
            .map(|(index, _)| index);
        let Some(first_uncensored_index) = first_uncensored else {
            return Ok(CensoredSdContext {
                observations: Vec::with_capacity(0),
                covariate_statistics: Vec::with_capacity(0),
                unique_covariates: Vec::with_capacity(0),
                thresholds: Vec::with_capacity(0),
            });
        };
        let capacity_upper_bound = n - first_uncensored_index;
        let mut obs = Vec::with_capacity(capacity_upper_bound);
        let mut thresholds = Vec::with_capacity(capacity_upper_bound);

        // Simultaneously deduplicate, copy over with index, and collect unique thresholds

        // First item
        let data_index = response_order[first_uncensored_index];
        thresholds.push(y[data_index]);
        obs.push(Observation {
            x: x[data_index],
            y: 0,
            observed: true,
            weight: 0.0, // placeholder; finalized from `last_w_accum`
        });
        debug_assert!(observed[data_index]);
        // Accumulator for the in-progress observation's weight in the caller's W precision.
        // Invariant: while the loop is running, `last_w_accum` holds the running W-precision
        // sum for `obs.last()`. Its `weight` field carries a 0.0 placeholder that is replaced
        // by the narrowed accumulator when the observation is finalized — either when a new
        // observation is pushed, or after the loop ends.
        let mut last_w_accum: W = weights[data_index];
        // Remaining items
        for &data_index in &response_order[first_uncensored_index + 1..] {
            let response_equal = y[data_index] == *thresholds.last().unwrap();
            let last_observation = obs.last().unwrap();
            let censoring_equal = observed[data_index] == last_observation.observed;
            let covariate_equal = x[data_index] == last_observation.x;

            let is_duplicate = response_equal && censoring_equal && covariate_equal;
            let is_following_censored =
                observed[data_index] && !last_observation.observed && covariate_equal;
            if is_duplicate || is_following_censored {
                last_w_accum = last_w_accum + weights[data_index];
            } else {
                // Finalize the previous in-progress observation: narrow once.
                obs.last_mut().unwrap().weight = last_w_accum.to_f32().unwrap();
                if !response_equal && observed[data_index] {
                    thresholds.push(y[data_index]);
                }
                obs.push(Observation {
                    x: x[data_index],
                    // If an observation is censored, we point to the previous (lower) threshold value here
                    y: thresholds.len() - 1,
                    observed: observed[data_index],
                    weight: 0.0, // placeholder; finalized from `last_w_accum`
                });
                last_w_accum = weights[data_index];
            }
        }
        // Finalize the last in-progress observation.
        obs.last_mut().unwrap().weight = last_w_accum.to_f32().unwrap();

        (obs, thresholds)
    };

    // Order of the observations (can sort unstable, we keep the same response order)
    let covariate_order = argsort_unstable_by::<Increasing, _>(
        |i, j| {
            observations_response_sorted[i]
                .x
                .total_cmp(&observations_response_sorted[j].x)
        },
        observations_response_sorted.len(),
    );

    let mut observations: Vec<Observation<usize, usize, bool, f32>> = vec![
        Observation {
            x: usize::MAX,
            y: usize::MAX,
            observed: false,
            weight: f32::NAN,
        };
        observations_response_sorted.len()
    ];
    let capacity_upper_bound = observations_response_sorted.len();
    let mut unique_covariates = Vec::with_capacity(capacity_upper_bound);
    let mut covariate_statistics = Vec::with_capacity(capacity_upper_bound);

    // First item
    let data_index = covariate_order[0];
    let observation = &observations_response_sorted[data_index];
    observations[data_index] = Observation {
        x: 0,
        y: observation.y,
        observed: observation.observed,
        weight: observation.weight,
    };
    unique_covariates.push(observation.x);
    covariate_statistics.push(CovariateStatistic {
        weight: observation.weight,
        cumulative_weight: 0.0,
    });
    // Remaining items
    for data_index in covariate_order.into_iter().skip(1) {
        let observation = &observations_response_sorted[data_index];
        if observation.x != *unique_covariates.last().unwrap() {
            unique_covariates.push(observation.x);
            let last_statistic = covariate_statistics.last_mut().unwrap();
            last_statistic.cumulative_weight += last_statistic.weight;
            covariate_statistics.push(CovariateStatistic {
                weight: 0.0,
                cumulative_weight: covariate_statistics.last().unwrap().cumulative_weight,
            });
        }
        observations[data_index] = Observation {
            x: unique_covariates.len() - 1,
            y: observation.y,
            observed: observation.observed,
            weight: observation.weight,
        };
        covariate_statistics.last_mut().unwrap().weight += observation.weight;
    }
    let last_statistic = covariate_statistics.last_mut().unwrap();
    last_statistic.cumulative_weight += last_statistic.weight;
    debug_assert!(observations.iter().all(|o| !o.weight.is_nan()));

    Ok(CensoredSdContext {
        observations,
        covariate_statistics,
        unique_covariates,
        thresholds,
    })
}

#[cfg(test)]
mod test {
    use crate::structures::Observation;
    use crate::total_order::stochastic_dominance::censored::preprocessing::preprocess;
    use crate::total_order::stochastic_dominance::censored::structures::CensoredSdContext;
    use crate::total_order::structures::CovariateStatistic;

    #[test]
    fn test_trivial_single_observation() {
        let context = preprocess(&[5.0], &[6.0], &[true], &[2.0]).ok().unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![Observation {
                    x: 0,
                    y: 0,
                    observed: true,
                    weight: 2.0,
                },],
                covariate_statistics: vec![CovariateStatistic {
                    weight: 2.0,
                    cumulative_weight: 2.0,
                },],
                unique_covariates: vec![5.0],
                thresholds: vec![6.0],
            }
        );
    }

    #[test]
    fn test_trivial_single_covariate() {
        let context = preprocess(&[5.0, 5.0], &[6.5, 6.0], &[true, true], &[1.0, 2.0])
            .ok()
            .unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 2.0,
                    },
                    Observation {
                        x: 0,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![CovariateStatistic {
                    weight: 3.0,
                    cumulative_weight: 3.0,
                },],
                unique_covariates: vec![5.0],
                thresholds: vec![6.0, 6.5],
            }
        );
    }

    #[test]
    fn test_trivial_single_response() {
        let context = preprocess(&[5.0, 7.0], &[6.0, 6.0], &[true, true], &[1.0, 3.0])
            .ok()
            .unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 1,
                        y: 0,
                        observed: true,
                        weight: 3.0,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 1.0,
                    },
                    CovariateStatistic {
                        weight: 3.0,
                        cumulative_weight: 4.0,
                    },
                ],
                unique_covariates: vec![5.0, 7.0],
                thresholds: vec![6.0],
            }
        );
    }

    #[test]
    fn test_monotone() {
        let context = preprocess(
            &[2.0, 1.0, 3.0],
            &[2.0, 1.0, 3.0],
            &[false, true, true],
            &[1.0, 9.0, 1.0],
        )
        .ok()
        .unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 9.0,
                    },
                    Observation {
                        x: 1,
                        y: 0,
                        observed: false,
                        weight: 1.0,
                    },
                    Observation {
                        x: 2,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![
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
                unique_covariates: vec![1.0, 2.0, 3.0],
                thresholds: vec![1.0, 3.0],
            }
        );
    }

    #[test]
    fn test_duplicate_response() {
        let context =
            preprocess(&[1.0, 3.0, 2.0], &[3.0, 3.0, 5.0], &[true; 3], &[1.0; 3]).unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 2,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 1,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 1.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 2.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 3.0,
                    },
                ],
                unique_covariates: vec![1.0, 2.0, 3.0],
                thresholds: vec![3.0, 5.0],
            },
        );
    }

    #[test]
    fn test_not_monotone_censored_2() {
        let context = preprocess(&[1.0, 2.0], &[4.0, 3.0], &[true, false], &[1.0, 1.0]).unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![Observation {
                    x: 0,
                    y: 0,
                    observed: true,
                    weight: 1.0,
                },],
                covariate_statistics: vec![CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 1.0,
                },],
                unique_covariates: vec![1.0],
                thresholds: vec![4.0],
            }
        );
    }

    #[test]
    fn test_not_monotone_3() {
        let context = preprocess(
            &[2.0, 1.0, 3.0],
            &[3.0, 1.0, 2.0],
            &[false, true, true],
            &[1.0, 9.0, 1.0],
        )
        .unwrap();
        assert_eq!(
            context,
            CensoredSdContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 9.0,
                    },
                    Observation {
                        x: 2,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 1,
                        y: 1,
                        observed: false,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![
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
                unique_covariates: vec![1.0, 2.0, 3.0],
                thresholds: vec![1.0, 2.0],
            }
        );
    }
}
