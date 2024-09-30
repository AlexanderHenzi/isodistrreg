use crate::error::Error;
use crate::progress::ProgressTracker;
use crate::routines::kaplan_meier;
use crate::structures::{Decreasing, Direction, Increasing};
use crate::total_order::hazard_rate_order::preprocessing::preprocess_censored;
use crate::total_order::routines::{pool_partitions_from_right, single_response};
use crate::total_order::structures::{AlgorithmContext, AlgorithmOutput, WeightedPartition};
use std::iter::{repeat, repeat_n};
use std::mem;

pub fn algorithm<D: Direction>(
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weight: &[f64],
    epsilon: f64,
    progress: &dyn ProgressTracker,
) -> Result<AlgorithmOutput, Error> {
    if !D::IS_INCREASING {
        unimplemented!();
    }

    let AlgorithmContext {
        observations,
        mut covariate_statistics,
        unique_responses,
        unique_covariates,
    } = preprocess_censored(x, y, observed, weight);
    progress.set_total(unique_responses.len());
    let n_covariate = covariate_statistics.len();

    if unique_responses.len() == 1 {
        // Single threshold -> simple binary isotonic regression with censoring amount
        return Ok(AlgorithmOutput {
            cdfs: single_response::<D, _>(observations, covariate_statistics),
            unique_covariates,
            thresholds: unique_responses,
        });
    }
    if covariate_statistics.len() == 1 {
        // Single covariate -> a single empirical cdf
        // Observations have been deduplicated so responses are unique
        return Ok(AlgorithmOutput {
            cdfs: kaplan_meier(observations.iter().copied(), covariate_statistics[0].weight),
            unique_covariates,
            thresholds: unique_responses,
        });
    }

    // At least two thresholds
    let mut cdfs = Vec::with_capacity(unique_responses.len() * n_covariate);
    let mut partitions = Vec::with_capacity(n_covariate);

    let mut data_index = 0;
    let mut threshold = observations[data_index].y;
    while data_index < observations.len() && !observations[data_index].observed {
        let observation = &observations[data_index];
        covariate_statistics[observation.x].weight -= observation.weight;
        data_index += 1;
        if data_index < observations.len() && observations[data_index].y != threshold {
            cdfs.extend(repeat_n(0.0, n_covariate));
            threshold = observations[data_index].y;
        }
    }
    if data_index == observations.len() {
        cdfs.extend(repeat_n(0.0, n_covariate));
        return Ok(AlgorithmOutput {
            cdfs,
            unique_covariates,
            thresholds: unique_responses,
        });
    }

    let mut estimators = vec![1.0; n_covariate];
    let mut at_risk: Vec<_> = covariate_statistics.iter().map(|cs| cs.weight).collect();

    let mut zero_count = 0;

    let mut survival: Vec<_> = {
        // Compute hazard rates in a sparse way
        let first_observation = &observations[data_index];
        assert!(first_observation.observed);
        estimators[first_observation.x] -= first_observation.weight / at_risk[first_observation.x];
        let total_at_risk = at_risk[..=first_observation.x].iter().sum();
        partitions.push(WeightedPartition {
            index: first_observation.x + 1,
            weight: total_at_risk,
            value: first_observation.weight / total_at_risk,
        });
        data_index += 1;
        covariate_statistics[first_observation.x].weight -= first_observation.weight;

        while observations[data_index].y == threshold && observations[data_index].observed {
            let previous = &observations[data_index - 1];
            let observation = &observations[data_index];

            estimators[observation.x] -= observation.weight / at_risk[observation.x];
            let total_at_risk = at_risk[(previous.x + 1)..=observation.x].iter().sum();
            partitions.push(WeightedPartition {
                index: observation.x + 1,
                weight: total_at_risk,
                value: observation.weight / total_at_risk,
            });
            pool_partitions_from_right::<Decreasing>(&mut partitions);

            covariate_statistics[observation.x].weight -= observation.weight;
            data_index += 1;
        }

        if partitions[0].value >= 1.0 - epsilon {
            zero_count = partitions[0].index;
        }

        let values = partitions
            .drain(..)
            .scan(0, |start_index, partition| {
                let previous_index = mem::replace(start_index, partition.index);
                Some(repeat_n(partition.value, partition.index - previous_index))
            })
            .flatten()
            .chain(repeat(0.0))
            .take(n_covariate);
        cdfs.extend(values);
        progress.increment();
        cdfs.iter().map(|v| 1.0 - v).collect()
    };
    for i in 0..n_covariate {
        at_risk[i] *= survival[i];
    }
    while observations[data_index].y == threshold {
        let observation = observations[data_index];
        debug_assert!(!observation.observed);

        at_risk[observation.x] *=
            1.0 - observation.weight / covariate_statistics[observation.x].weight;
        covariate_statistics[observation.x].weight -= observation.weight;
        data_index += 1;
    }

    while data_index < observations.len() {
        threshold = observations[data_index].y;

        // Update zero count
        while data_index < observations.len() && observations[data_index].y == threshold {
            while survival[zero_count] <= 0.0 + epsilon {
                zero_count += 1;
            }

            let observation = &observations[data_index];
            // TODO: Numerics
            if observation.x == zero_count
                && observation.observed
                && (observation.weight - covariate_statistics[observation.x].weight).abs()
                    <= epsilon
            {
                // No need to update the weight total in the covariate static or estimator, they
                // won't be used from now on
                zero_count += 1;
                // Skip the estimators already zero
                while zero_count < n_covariate && estimators[zero_count] <= 0.0 + epsilon {
                    zero_count += 1;
                }

                data_index += 1;
            } else {
                break;
            }
        }

        // Update non-zero items
        let mut i = zero_count;
        while data_index < observations.len()
            && observations[data_index].y == threshold
            && observations[data_index].observed
        {
            let observation = &observations[data_index];

            let share_of_remaining =
                observation.weight / covariate_statistics[observation.x].weight;
            estimators[observation.x] *= 1.0 - share_of_remaining;
            covariate_statistics[observation.x].weight -= observation.weight;

            partitions.push(WeightedPartition {
                index: i + 1,
                weight: at_risk[i],
                value: estimators[i] / survival[i],
            });
            i += 1;
            while i < observation.x {
                partitions.push(WeightedPartition {
                    index: i + 1,
                    weight: at_risk[i],
                    value: estimators[i] / survival[i],
                });
                pool_partitions_from_right::<Increasing>(&mut partitions);
                i += 1;
            }
            partitions.push(WeightedPartition {
                index: observation.x + 1,
                weight: at_risk[observation.x],
                value: estimators[observation.x] / survival[observation.x],
            });
            pool_partitions_from_right::<Increasing>(&mut partitions);
            i += 1;

            data_index += 1;
        }
        while i < n_covariate {
            partitions.push(WeightedPartition {
                index: i + 1,
                weight: at_risk[i],
                value: estimators[i] / survival[i],
            });
            pool_partitions_from_right::<Increasing>(&mut partitions);
            i += 1;
        }
        while data_index < observations.len() && observations[data_index].y == threshold {
            let observation = observations[data_index];
            debug_assert!(!observation.observed);

            at_risk[observation.x] *=
                1.0 - observation.weight / covariate_statistics[observation.x].weight;
            covariate_statistics[observation.x].weight -= observation.weight;
            data_index += 1;
        }

        // Save results
        cdfs.extend(repeat_n(1.0, zero_count));
        let mut start_index = zero_count;
        for partition in partitions.drain(..) {
            for index in start_index..partition.index {
                // Previous iteration update
                at_risk[index] *= partition.value;
                // Current iteration update
                survival[index] *= partition.value;
                // Write out result, clamp to ensure numerical noise doesn't get us out of [0, 1]
                cdfs.push(1.0 - survival[index].clamp(0.0, 1.0));
            }
            start_index = partition.index;
        }
        progress.increment();
    }

    Ok(AlgorithmOutput {
        cdfs,
        unique_covariates,
        thresholds: unique_responses,
    })
}

#[cfg(test)]
mod test {
    use crate::structures::Increasing;
    use crate::test::is_relative_eq_vec;
    use crate::total_order::hazard_rate_order::censored::algorithm;
    use crate::total_order::structures::AlgorithmOutput;

    fn execute_test<const N: usize, const N_COVARIATE: usize, const N_THRESHOLD: usize>(
        x: [f64; N],
        y: [f64; N],
        observed: [bool; N],
        weight: [f64; N],
        expected: [[f64; N_COVARIATE]; N_THRESHOLD],
    ) {
        let AlgorithmOutput {
            cdfs,
            unique_covariates,
            thresholds,
        } = algorithm::<Increasing>(&x, &y, &observed, &weight, 1e-10, &crate::NoProgress)
            .ok()
            .unwrap();
        assert_eq!(unique_covariates.len(), N_COVARIATE);
        assert_eq!(thresholds.len(), N_THRESHOLD);
        let expected_flat: Vec<_> = expected.iter().flatten().copied().collect();
        assert!(
            is_relative_eq_vec(&cdfs, &expected_flat),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected_flat,
        );
    }

    #[test]
    fn test_weighted_duplicate() {
        execute_test(
            [0.0, 2.0, 1.0, 1.0, 3.0],
            [6.0, 7.0, 8.0, 8.0, 9.0],
            [true; 5],
            [1.0, 2.0, 1.0, 1.0, 1.0],
            [
                [1.0, 0.0, 0.0, 0.0],
                [1.0, 0.5, 0.5, 0.0],
                [1.0, 1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0, 1.0],
            ],
        );
    }

    #[test]
    fn test_mixed_tonicity_4() {
        execute_test(
            [0.0, 1.0, 2.0, 3.0],
            [2.0, 1.0, 3.0, 4.0],
            [true; 4],
            [1.0; 4],
            [
                [0.5, 0.5, 0.0, 0.0],
                [1.0, 1.0, 0.0, 0.0],
                [1.0, 1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0, 1.0],
            ],
        );
    }

    #[test]
    fn test_mixed_tonicity_3() {
        execute_test(
            [0.0, 1.0, 2.0],
            [3.0, 1.0, 2.0],
            [true; 3],
            [1.0; 3],
            [[0.5, 0.5, 0.0], [0.75, 0.75, 0.5], [1.0, 1.0, 1.0]],
        );
    }

    #[test]
    fn test_mixed_tonicity_6() {
        execute_test(
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            [6.0, 2.0, 3.0, 4.0, 1.0, 5.0],
            [true; 6],
            [1.0; 6],
            [
                [1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 0.0],
                [1.0 / 2.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
                [7.0 / 10.0, 7.0 / 10.0, 3.0 / 5.0, 1.0 / 2.0, 1.0 / 2.0, 0.0],
                [
                    17.0 / 20.0,
                    17.0 / 20.0,
                    4.0 / 5.0,
                    3.0 / 4.0,
                    3.0 / 4.0,
                    0.0,
                ],
                [
                    37.0 / 40.0,
                    37.0 / 40.0,
                    9.0 / 10.0,
                    7.0 / 8.0,
                    7.0 / 8.0,
                    0.5,
                ],
                [1.0; 6],
            ],
        );
    }

    #[test]
    fn test_5() {
        execute_test(
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [3.0, 2.0, 1.0, 5.0, 4.0],
            [true, false, false, false, false],
            [1.0; 5],
            [[1.0, 0.0, 0.0, 0.0, 0.0]],
        );
    }

    #[test]
    fn test_5_1() {
        execute_test(
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [3.0, 2.0, 1.0, 5.0, 4.0],
            [false, false, true, false, false],
            [1.0; 5],
            [[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0]],
        );
    }
}
