use crate::routines::argsort_unstable_by;
use crate::structures::{Increasing, Observation};
use crate::total_order::structures::AlgorithmContext;
use crate::total_order::structures::CovariateStatistic;

pub fn preprocess_uncensored(x: &[f64], y: &[f64], weight: &[f64]) -> AlgorithmContext<()> {
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

pub fn preprocess<S: Ord, F: Fn(usize) -> S>(
    x: &[f64],
    y: &[f64],
    observed: F,
    weight: &[f64],
) -> AlgorithmContext<S> {
    let n = x.len();

    // Determine the order by covariate and response
    let order = argsort_unstable_by::<Increasing, _>(
        |i, j| {
            x[i].total_cmp(&x[j])
                .then(y[i].total_cmp(&y[j]))
                .then(observed(i).cmp(&observed(j)).reverse())
        },
        n,
    );

    // While copying the data in the sorted order and aggregating identical observations by weight,
    // we track the unique covariates we encounter.
    let mut observations = Vec::with_capacity(n);
    let mut unique_covariates = Vec::with_capacity(n);
    let mut covariate_statistics = Vec::with_capacity(n);

    let mut covariate_sorted_indices = order.into_iter();
    {
        let index = covariate_sorted_indices.next().unwrap();
        observations.push(Observation {
            x: 0,
            y: y[index],
            observed: observed(index),
            weight: weight[index],
        });
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
            // At the same observation -> just accumulate observation weight
            last_observation.weight += weight[index];
        } else if covariate_equal {
            // New observation but same covariate -> accumulate covariate weight from observation
            // and add new observation
            current_statistic.weight += last_observation.weight;
            observations.push(Observation {
                x: unique_covariates.len() - 1, // we stay with the same covariate
                y: y[index],
                observed: observed(index),
                weight: weight[index],
            });
        } else {
            // New observation and new covariate -> accumulate covariate weight from observation,
            // add a new observation and collect the cumulative covariate weight and add a new
            // covariate statistic
            current_statistic.weight += last_observation.weight;
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
                weight: weight[index],
            });
            unique_covariates.push(x[index]);

            current_statistic = covariate_statistics.last_mut().unwrap();
        }
    }
    current_statistic.weight += observations.last().unwrap().weight;
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
