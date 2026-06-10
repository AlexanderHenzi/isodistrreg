use crate::progress::ProgressTracker;
use crate::routines::empirical_cdf;
use crate::structures::{Decreasing, Direction, Increasing};
use crate::total_order::routines::pool_partitions_from_right;
use crate::total_order::structures::{AlgorithmContext, WeightedPartition};
use crate::total_order::weight_noise_floor;
use std::iter::{repeat, repeat_n};
use std::mem;

/// Compute the IDR for a totally ordered covariate with uncensored hazard rate ordered responses.
///
/// # Arguments
///
/// - `context`: pre-built algorithm context (preprocessed input)
///
/// # Returns
///
/// A threshold-major matrix of cdf values.
///
/// # Complexity
///
/// TODO
///
/// # Examples
///
/// ```rust
/// use isodistrreg::{Increasing, NoProgress};
/// use isodistrreg::total_order::hazard_rate_order::{preprocess_uncensored, uncensored};
///
/// let covariates = [1.0, 2.0, 3.0];
/// let responses = [8.0, 6.0, 7.0];
/// let context = preprocess_uncensored(&covariates, &responses, &[1.0; 3]);
/// let cdfs = uncensored::<Increasing, _, _>(&context, &NoProgress);
///
/// let expected = vec![
///     0.5, 0.5, 0.0,
///     5.0 / 6.0, 5.0 / 6.0, 2.0 / 3.0,
///     1.0, 1.0, 1.0,
/// ];
/// for (r, e) in cdfs.into_iter().zip(expected.into_iter()) {
///     assert!((r - e).abs() < 1e-6);
/// }
/// ```
pub fn algorithm<D: Direction, X: crate::Float, Y: crate::Float>(
    context: &AlgorithmContext<X, Y, ()>,
    progress: &dyn ProgressTracker,
) -> Vec<f32> {
    if !D::IS_INCREASING {
        unimplemented!();
    }

    let AlgorithmContext {
        observations,
        covariate_statistics,
        unique_responses,
        unique_covariates,
    } = context;
    progress.set_total(unique_responses.len());
    let n_covariate = covariate_statistics.len();
    // Dynamic f32-noise floor for the running weight subtractions and survival/CDF saturation
    // checks below. The hazard-rate uncensored kernel has a slightly different "noise source"
    // from the SD-censored case (here we sum `survival[i] * survival[i] *
    // covariate_statistics[i].weight` and compare `estimators[i]` against 0/1) but the
    // Wilkinson-style bound in `weight_noise_floor` is conservative enough to cover both.
    let epsilon = weight_noise_floor(covariate_statistics.iter().map(|s| s.weight).sum());

    if unique_responses.len() == 1 {
        // Single threshold -> all one's
        return vec![1.0; n_covariate];
    }
    if unique_covariates.len() == 1 {
        // Single covariate -> a single empirical cdf
        // Observations have been deduplicated so responses are unique
        return empirical_cdf(observations.iter().copied(), covariate_statistics[0].weight);
    }

    // At least two thresholds
    let mut cdfs = Vec::with_capacity(unique_responses.len() * n_covariate);
    let mut partitions: Vec<WeightedPartition> = Vec::with_capacity(n_covariate);

    let mut estimators = vec![1.0; n_covariate];

    let mut data_index = 0;
    let mut zero_count = 0;

    let mut threshold = observations[data_index].y;

    let mut survival: Vec<_> = {
        // Compute hazard rates in a sparse way
        let first_observation = &observations[data_index];
        estimators[first_observation.x] -=
            first_observation.weight / covariate_statistics[first_observation.x].weight;
        let total_weight = covariate_statistics[..=first_observation.x]
            .iter()
            .map(|cs| cs.weight)
            .sum();
        partitions.push(WeightedPartition {
            index: first_observation.x + 1,
            weight: total_weight,
            value: first_observation.weight / total_weight,
        });
        data_index += 1;

        while observations[data_index].y == threshold {
            let previous = &observations[data_index - 1];
            let observation = &observations[data_index];

            estimators[observation.x] -=
                observation.weight / covariate_statistics[observation.x].weight;
            let total_weight = covariate_statistics[(previous.x + 1)..=observation.x]
                .iter()
                .map(|cs| cs.weight)
                .sum();
            partitions.push(WeightedPartition {
                index: observation.x + 1,
                weight: total_weight,
                value: observation.weight / total_weight,
            });
            pool_partitions_from_right::<Decreasing, _>(&mut partitions);

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
        // Write out results as a CDF (we computed survival), clamp to ensure numerical noise
        // doesn't get us out of [0, 1]
        cdfs.iter().map(|v| 1.0 - v.clamp(0.0, 1.0)).collect()
    };

    let last_threshold = observations.last().unwrap().y;
    while observations[data_index].y < last_threshold {
        threshold = observations[data_index].y;

        // Update zero count
        while observations[data_index].y == threshold {
            let observation = &observations[data_index];
            // TODO: Numerics
            if observation.x == zero_count
                && (observation.weight
                    - estimators[observation.x] * covariate_statistics[observation.x].weight)
                    .abs()
                    <= epsilon
            {
                // No need to update the weight total in the covariate static or estimator, they
                // won't be used from now on
                zero_count += 1;
                // Skip the estimators already zero
                while estimators[zero_count] <= 0.0 + epsilon {
                    zero_count += 1;
                }

                data_index += 1;
            } else {
                break;
            }
        }

        // Update non-zero items
        let mut i = zero_count;
        while observations[data_index].y == threshold {
            let observation = &observations[data_index];

            estimators[observation.x] -=
                observation.weight / covariate_statistics[observation.x].weight;

            while i <= observation.x {
                partitions.push(WeightedPartition {
                    index: i + 1,
                    weight: survival[i] * survival[i] * covariate_statistics[i].weight,
                    value: estimators[i] / survival[i],
                });
                pool_partitions_from_right::<Increasing, _>(&mut partitions);
                i += 1;
            }

            data_index += 1;
        }
        while i < n_covariate {
            partitions.push(WeightedPartition {
                index: i + 1,
                weight: survival[i] * survival[i] * covariate_statistics[i].weight,
                value: estimators[i] / survival[i],
            });
            pool_partitions_from_right::<Increasing, _>(&mut partitions);
            i += 1;
        }

        // Save results
        cdfs.extend(repeat_n(1.0, zero_count));
        let mut start_index = zero_count;
        for partition in partitions.drain(..) {
            for s in &mut survival[start_index..partition.index] {
                *s *= partition.value;
                cdfs.push(1.0 - *s);
            }
            start_index = partition.index;
        }
        progress.increment();
    }

    // Last iteration is all ones
    cdfs.extend(repeat_n(1.0, n_covariate));
    progress.increment();

    cdfs
}

#[cfg(test)]
mod test {
    use crate::preprocessing::validate;
    use crate::structures::Increasing;
    use crate::test::is_relative_eq_vec;
    use crate::total_order::hazard_rate_order::uncensored::algorithm;
    use crate::total_order::preprocessing::preprocess_uncensored;

    fn execute_test<const N: usize, const N_COVARIATE: usize, const N_THRESHOLD: usize>(
        x: [f64; N],
        y: [f64; N],
        weight: [f64; N],
        expected: [[f64; N_COVARIATE]; N_THRESHOLD],
    ) {
        let expected_flat: Vec<_> = expected.iter().flatten().copied().collect();
        validate(x.chunks_exact(1), &y, None, Some(&weight)).unwrap();
        let context = preprocess_uncensored(&x, &y, &weight);
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.unique_responses.len(), N_THRESHOLD);
        // The hazard-rate algorithm now runs in f32; narrow the test's f64 expected.
        let expected_f32: Vec<_> = expected_flat.iter().map(|&v| v as f32).collect();
        assert!(
            is_relative_eq_vec(&cdfs, &expected_f32),
            "Result:   {:?}\nExpected: {:?}\n",
            cdfs,
            expected,
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
                [1.0, 1.0, 1.0, 1.0],
            ],
        );
    }

    #[test]
    fn test_mixed_tonicity_4() {
        execute_test(
            [0.0, 1.0, 2.0, 3.0],
            [2.0, 1.0, 3.0, 4.0],
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
            [1.0; 3],
            [
                [0.5, 0.5, 0.0],
                [5.0 / 6.0, 5.0 / 6.0, 2.0 / 3.0],
                [1.0, 1.0, 1.0],
            ],
        );
    }

    #[test]
    fn test_mixed_tonicity_6() {
        execute_test(
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            [6.0, 2.0, 3.0, 4.0, 1.0, 5.0],
            [1.0; 6],
            [
                [1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0, 0.0],
                [1.0 / 2.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0],
                [
                    25.0 / 34.0,
                    25.0 / 34.0,
                    11.0 / 17.0,
                    1.0 / 2.0,
                    1.0 / 2.0,
                    0.0,
                ],
                // TODO: Manually verify these test values
                [
                    803.0 / 884.0,
                    803.0 / 884.0,
                    194.0 / 221.0,
                    43.0 / 52.0,
                    43.0 / 52.0,
                    0.0,
                ],
                [
                    846499.0 / 853060.0,
                    846499.0 / 853060.0,
                    211078.0 / 213265.0,
                    49451.0 / 50180.0,
                    49451.0 / 50180.0,
                    884.0 / 965.0,
                ],
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            ],
        );
    }
}
