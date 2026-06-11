use crate::progress::ProgressTracker;
use crate::routines::kaplan_meier;
use crate::structures::{Decreasing, Direction, Increasing};
use crate::total_order::routines::{pool_partitions_from_right, single_response};
use crate::total_order::structures::{AlgorithmContext, WeightedPartition};
use itertools::Itertools;
use std::iter::{repeat, repeat_n};
use std::mem;

pub fn algorithm<D: Direction, X: crate::Float, Y: crate::Float>(
    context: &AlgorithmContext<X, Y, bool>,
    progress: &dyn ProgressTracker,
) -> Vec<f32> {
    if !D::IS_INCREASING {
        unimplemented!();
    }

    let observations = &context.observations;
    // The algorithm mutates `covariate_statistics` (decrementing weights as observations are
    // consumed). Clone once at entry so the caller's context stays immutable.
    let mut covariate_statistics = context.covariate_statistics.clone();
    let unique_responses = &context.unique_responses;
    progress.set_total(unique_responses.len());
    let n_covariate = covariate_statistics.len();

    if unique_responses.len() == 1 {
        // Single threshold -> simple binary isotonic regression with censoring amount
        return single_response::<D, _>(observations.clone(), covariate_statistics);
    }
    if covariate_statistics.len() == 1 {
        // Single covariate -> a single Kaplan-Meier curve. `kaplan_meier` emits one
        // value per distinct response INCLUDING censored-only ones, but the fit's
        // thresholds (`unique_responses`) keep only responses with an observed event —
        // restrict the curve to those.
        let km = kaplan_meier(observations.iter().copied(), covariate_statistics[0].weight);
        return observations
            .iter()
            .map(|o| o.y)
            .dedup()
            .zip(km)
            .filter(|(y, _)| {
                unique_responses
                    .binary_search_by(|t| t.total_cmp(y))
                    .is_ok()
            })
            .map(|(_, v)| v)
            .collect();
    }

    // At least two thresholds
    let mut cdfs = Vec::with_capacity(unique_responses.len() * n_covariate);
    let mut partitions: Vec<WeightedPartition> = Vec::with_capacity(n_covariate);

    // Index (+1; 0 = none) of each covariate's last positive-weight observation. Once
    // it is consumed no mass remains in the group, so the raw Kaplan-Meier estimator
    // (when it is an event) respectively the at-risk mass (when censored) is exactly
    // 0 — a purely combinatorial fact. Pin those cases below: the running f32 weight
    // bookkeeping drifts a few ulps off the exact zero otherwise. Whether the last
    // observation is an event decides DEATH (CDF pinned at 1) versus FROZEN (all mass
    // censored away, estimator stays positive, CDF stays below 1).
    let mut last_positive = vec![0usize; n_covariate];
    for (d, o) in observations.iter().enumerate() {
        if o.weight > 0.0 {
            last_positive[o.x] = d + 1;
        }
    }

    let mut data_index = 0;
    // Leading censored observations (before the first event) only shrink the
    // available mass; their response values are not thresholds (no event), so no CDF
    // row is emitted for them.
    while data_index < observations.len() && !observations[data_index].observed {
        let observation = &observations[data_index];
        covariate_statistics[observation.x].weight -= observation.weight;
        if data_index + 1 == last_positive[observation.x] {
            // Fully censored group: no positive mass remains — exactly 0, so the
            // at-risk vector below starts at an exact 0 instead of subtraction drift.
            covariate_statistics[observation.x].weight = 0.0;
        }
        data_index += 1;
    }
    if data_index == observations.len() {
        // No events at all: there are no thresholds and the fit is the empty sub-CDF.
        debug_assert!(unique_responses.is_empty());
        return cdfs;
    }
    let mut threshold = observations[data_index].y;

    let mut estimators = vec![1.0; n_covariate];
    let mut at_risk: Vec<_> = covariate_statistics.iter().map(|cs| cs.weight).collect();

    let mut zero_count = 0;

    let mut survival: Vec<_> = {
        // Compute hazard rates in a sparse way
        let first_observation = &observations[data_index];
        assert!(first_observation.observed);
        // Zero-weight events carry no mass (and would divide 0/0 against an exactly
        // pinned at-risk mass).
        if first_observation.weight > 0.0 {
            estimators[first_observation.x] -=
                first_observation.weight / at_risk[first_observation.x];
        }
        if data_index + 1 == last_positive[first_observation.x] {
            // This event consumes the group's whole remaining mass — death, exactly.
            estimators[first_observation.x] = 0.0;
        }
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

            if observation.weight > 0.0 {
                estimators[observation.x] -= observation.weight / at_risk[observation.x];
            }
            if data_index + 1 == last_positive[observation.x] {
                // This event consumes the group's whole remaining mass — death, exactly.
                estimators[observation.x] = 0.0;
            }
            let total_at_risk = at_risk[(previous.x + 1)..=observation.x].iter().sum();
            partitions.push(WeightedPartition {
                index: observation.x + 1,
                weight: total_at_risk,
                value: observation.weight / total_at_risk,
            });
            pool_partitions_from_right::<Decreasing, _>(&mut partitions);

            covariate_statistics[observation.x].weight -= observation.weight;
            data_index += 1;
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
        cdfs.iter().map(|v| 1.0 - v.clamp(0.0, 1.0)).collect()
    };
    // Pin every leading group that can no longer hold mass — exactly, no epsilon:
    // - `estimators == 0.0` marks groups whose last positive observation was an event
    //   (snapped above via `last_positive`): they died at the first threshold.
    // - `survival == 0.0` additionally catches fully-censored covariates lying INSIDE
    //   a dying event's partition span: their at-risk mass is exactly 0, so the span's
    //   value is exactly 1 and their fitted survival exactly 0. Without pinning, the
    //   next threshold would compute the ratio `estimators / survival = 1/0`.
    // Groups whose last observation is censored and that sit outside any dying span
    // keep a positive estimator and survival: they are frozen, not dead, and their
    // CDF must stay below 1 — neither check ever fires for them.
    while zero_count < n_covariate
        && (estimators[zero_count] == 0.0 || survival[zero_count] == 0.0)
    {
        zero_count += 1;
    }
    for i in 0..n_covariate {
        at_risk[i] *= survival[i];
    }
    while data_index < observations.len() && observations[data_index].y == threshold {
        let observation = observations[data_index];
        debug_assert!(!observation.observed);

        // Zero-weight observations carry no mass (and could divide 0/0 when the
        // covariate's remaining mass is already exhausted).
        if observation.weight > 0.0 {
            at_risk[observation.x] *=
                1.0 - observation.weight / covariate_statistics[observation.x].weight;
            covariate_statistics[observation.x].weight -= observation.weight;
            if data_index + 1 == last_positive[observation.x] {
                at_risk[observation.x] = 0.0;
            }
        }
        data_index += 1;
    }

    while data_index < observations.len() {
        threshold = observations[data_index].y;
        // Events sort before censored observations within a tied response, so this
        // threshold has at least one event iff its first observation is observed.
        // Event-free thresholds are not part of `unique_responses` and get no CDF
        // row; their censored observations only shrink the at-risk masses below.
        let emit_row = observations[data_index].observed;

        if emit_row {
            // Update zero count
            while data_index < observations.len() && observations[data_index].y == threshold {
                // Fitted survival reaches exactly 0 only for groups that died (their
                // snapped 0-estimator makes the partition value exactly 0) or whose
                // censored-away mass sat inside a dying span; frozen groups keep a
                // strictly positive survival and are never pinned.
                while zero_count < n_covariate && survival[zero_count] == 0.0 {
                    zero_count += 1;
                }

                let observation = &observations[data_index];
                // The group at the cursor dies here exactly when this event is its
                // last positive-weight observation: censored mass at earlier responses
                // already left, none remains at later ones, so the event consumes the
                // whole remaining mass — an exact index check, no weight comparison.
                if observation.x == zero_count
                    && observation.observed
                    && data_index + 1 == last_positive[observation.x]
                {
                    // No need to update the weight total in the covariate static or estimator,
                    // they won't be used from now on
                    zero_count += 1;
                    // Skip groups that already died earlier (estimators snapped to an
                    // exact 0 at their final event; frozen groups stay positive)
                    while zero_count < n_covariate && estimators[zero_count] == 0.0 {
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

                // Zero-weight events carry no mass (and would divide 0/0 once the
                // group's remaining mass is exactly pinned to 0).
                if observation.weight > 0.0 {
                    let share_of_remaining =
                        observation.weight / covariate_statistics[observation.x].weight;
                    estimators[observation.x] *= 1.0 - share_of_remaining;
                    covariate_statistics[observation.x].weight -= observation.weight;
                }
                if data_index + 1 == last_positive[observation.x] {
                    // This event consumes the group's whole remaining mass — death,
                    // exactly.
                    estimators[observation.x] = 0.0;
                }

                // Push every covariate up to and including the event's, pooling after
                // each push (mirrors the uncensored kernel). Pushing the cursor and the
                // event covariate separately would double-push when the event sits at
                // the cursor, skipping the next covariate and corrupting the row.
                while i <= observation.x {
                    partitions.push(WeightedPartition {
                        index: i + 1,
                        weight: at_risk[i],
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
                    weight: at_risk[i],
                    value: estimators[i] / survival[i],
                });
                pool_partitions_from_right::<Increasing, _>(&mut partitions);
                i += 1;
            }
        }
        while data_index < observations.len() && observations[data_index].y == threshold {
            let observation = observations[data_index];
            debug_assert!(!observation.observed);

            // Zero-weight observations carry no mass (and could divide 0/0 when the
            // covariate's remaining mass is already exhausted).
            if observation.weight > 0.0 {
                at_risk[observation.x] *=
                    1.0 - observation.weight / covariate_statistics[observation.x].weight;
                covariate_statistics[observation.x].weight -= observation.weight;
                if data_index + 1 == last_positive[observation.x] {
                    at_risk[observation.x] = 0.0;
                }
            }
            data_index += 1;
        }

        if emit_row {
            // Save results
            cdfs.extend(repeat_n(1.0, zero_count));
            let mut start_index = zero_count;
            for partition in partitions.drain(..) {
                for index in start_index..partition.index {
                    // Previous iteration update
                    at_risk[index] *= partition.value;
                    // Current iteration update
                    survival[index] *= partition.value;
                    // Write out result, clamp to ensure numerical noise doesn't get us out of
                    // [0, 1]
                    cdfs.push(1.0 - survival[index].clamp(0.0, 1.0));
                }
                start_index = partition.index;
            }
            progress.increment();
        }
    }

    cdfs
}

#[cfg(test)]
mod test {
    use crate::structures::Increasing;
    use crate::test::is_relative_eq_vec;
    use crate::total_order::hazard_rate_order::censored::algorithm;
    use crate::total_order::hazard_rate_order::preprocessing::preprocess_censored;

    fn execute_test<const N: usize, const N_COVARIATE: usize, const N_THRESHOLD: usize>(
        x: [f64; N],
        y: [f64; N],
        observed: [bool; N],
        weight: [f64; N],
        expected: [[f64; N_COVARIATE]; N_THRESHOLD],
    ) {
        let context = preprocess_censored(&x, &y, &observed, &weight);
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.unique_responses.len(), N_THRESHOLD);
        let expected_flat: Vec<_> = expected.iter().flatten().copied().collect();
        // The hazard-rate algorithm now runs in f32; narrow the test's f64 expected.
        let expected_f32: Vec<_> = expected_flat.iter().map(|&v| v as f32).collect();
        assert!(
            is_relative_eq_vec(&cdfs, &expected_f32),
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
