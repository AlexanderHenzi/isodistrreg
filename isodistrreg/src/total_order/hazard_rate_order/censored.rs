use crate::progress::ProgressTracker;
use crate::routines::kaplan_meier;
use crate::structures::{Decreasing, Direction, Increasing};
use crate::total_order::routines::{
    pool_partitions_from_right_zero_weight_blocks, single_response,
};
use crate::total_order::structures::{CensoredContext, WeightedPartition};
use std::iter::{repeat, repeat_n};
use std::mem;

pub fn algorithm<D: Direction, X: crate::Float, Y: crate::Float>(
    context: &CensoredContext<X, Y>,
    progress: &dyn ProgressTracker,
) -> Vec<f32> {
    if !D::IS_INCREASING {
        unimplemented!();
    }

    let observations = &context.observations;
    // The algorithm mutates `covariate_statistics` (decrementing weights as observations are
    // consumed). Clone once at entry so the caller's context stays immutable.
    let mut covariate_statistics = context.covariate_statistics.clone();
    progress.set_total(context.n_threshold());
    let n_covariate = covariate_statistics.len();

    if context.n_threshold() == 0 {
        // No events at all: there are no thresholds and the fit is the empty sub-CDF.
        debug_assert!(context.unique_covariates.is_empty());
        return Vec::with_capacity(0);
    }
    if context.n_threshold() == 1 {
        // Single threshold -> simple binary isotonic regression with censoring amount
        return single_response::<D, _>(observations.clone(), covariate_statistics);
    }
    if n_covariate == 1 {
        // Single covariate -> a single Kaplan-Meier curve. Censored responses are
        // snapped onto event thresholds during preprocessing, so the curve has exactly
        // one value per threshold.
        return kaplan_meier(observations.iter().copied(), covariate_statistics[0].weight);
    }

    // At least two thresholds
    let mut cdfs = Vec::with_capacity(context.n_threshold() * n_covariate);
    let mut partitions: Vec<WeightedPartition> = Vec::with_capacity(n_covariate);

    // Index (+1; 0 = none) of each covariate's last observation — preprocessing
    // guarantees every observation in the context has positive weight. Once it is
    // consumed no mass remains in the group, so the raw Kaplan-Meier estimator (when
    // it is an event) respectively the at-risk mass (when censored) is exactly 0 — a
    // purely combinatorial fact. Pin those cases below: the running f32 weight
    // bookkeeping drifts a few ulps off the exact zero otherwise. Whether the last
    // observation is an event decides DEATH (CDF pinned at 1) versus FROZEN (all mass
    // censored away, estimator stays positive, CDF stays below 1).
    let mut last_positive = vec![0usize; n_covariate];
    for (d, o) in observations.iter().enumerate() {
        debug_assert!(o.weight > 0.0);
        last_positive[o.x] = d + 1;
    }

    // Preprocessing discards censored observations sorting before the first event (they
    // only shrink mass that the covariate statistics already exclude), so the stream
    // starts with the first threshold's events.
    let mut data_index = 0;
    let mut threshold = observations[data_index].y;
    debug_assert!(observations[data_index].observed);

    let mut estimators = vec![1.0; n_covariate];
    let mut at_risk: Vec<_> = covariate_statistics.iter().map(|cs| cs.weight).collect();

    let mut zero_count = 0;

    let mut survival: Vec<_> = {
        // Compute hazard rates in a sparse way
        let first_observation = &observations[data_index];
        assert!(first_observation.observed);
        estimators[first_observation.x] -= first_observation.weight / at_risk[first_observation.x];
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

            estimators[observation.x] -= observation.weight / at_risk[observation.x];
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
            pool_partitions_from_right_zero_weight_blocks::<Decreasing, _>(&mut partitions);

            covariate_statistics[observation.x].weight -= observation.weight;
            data_index += 1;
        }

        let values = partitions
            .drain(..)
            .scan(0, |start_index, partition| {
                let previous_index = mem::replace(start_index, partition.index);
                // The raw first-threshold hazard can drift above 1 (the at-risk mass is
                // an f32 subtraction while the event weight is accumulated directly, e.g.
                // 0.2 / (10.2 - 10.0) = 1.000001); clamp into [0, 1] before repeating.
                // This is the earliest safe point: clamping at push time would feed
                // altered values into the pooling averages.
                Some(repeat_n(
                    partition.value.clamp(0.0, 1.0),
                    partition.index - previous_index,
                ))
            })
            .flatten()
            .chain(repeat(0.0))
            .take(n_covariate);
        cdfs.extend(values);
        progress.increment();
        cdfs.iter().map(|v| 1.0 - v).collect()
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
    while zero_count < n_covariate && (estimators[zero_count] == 0.0 || survival[zero_count] == 0.0)
    {
        zero_count += 1;
    }
    for i in 0..n_covariate {
        at_risk[i] *= survival[i];
    }
    while data_index < observations.len() && observations[data_index].y == threshold {
        let observation = observations[data_index];
        debug_assert!(!observation.observed);

        at_risk[observation.x] *=
            1.0 - observation.weight / covariate_statistics[observation.x].weight;
        covariate_statistics[observation.x].weight -= observation.weight;
        if data_index + 1 == last_positive[observation.x] {
            // Fully censored away: the remaining at-risk mass is exactly 0.
            at_risk[observation.x] = 0.0;
        }
        data_index += 1;
    }

    while data_index < observations.len() {
        threshold = observations[data_index].y;
        // Censored responses are snapped onto event thresholds during preprocessing
        // (events sorting before censored within a threshold), so every remaining
        // threshold starts with an event and emits a CDF row.
        debug_assert!(observations[data_index].observed);

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

            let share_of_remaining =
                observation.weight / covariate_statistics[observation.x].weight;
            estimators[observation.x] *= 1.0 - share_of_remaining;
            covariate_statistics[observation.x].weight -= observation.weight;
            if data_index + 1 == last_positive[observation.x] {
                // This event consumes the group's whole remaining mass — death,
                // exactly.
                estimators[observation.x] = 0.0;
            }

            // Push every covariate up to and including the event's, pooling after
            // each push (mirrors the uncensored kernel). Pushing the cursor and the
            // event covariate separately would double-push when the event sits at
            // the cursor, skipping the next covariate and corrupting the row.
            // The least-squares weight of the ratio regression is w·S(t-)²
            // (mirroring the uncensored kernel); `at_risk` already carries one
            // factor of the fitted survival, so multiply in the second.
            while i <= observation.x {
                partitions.push(WeightedPartition {
                    index: i + 1,
                    weight: at_risk[i] * survival[i],
                    value: estimators[i] / survival[i],
                });
                pool_partitions_from_right_zero_weight_blocks::<Increasing, _>(&mut partitions);
                i += 1;
            }

            data_index += 1;
        }
        while i < n_covariate {
            partitions.push(WeightedPartition {
                index: i + 1,
                weight: at_risk[i] * survival[i],
                value: estimators[i] / survival[i],
            });
            pool_partitions_from_right_zero_weight_blocks::<Increasing, _>(&mut partitions);
            i += 1;
        }

        // Censored observations at this threshold leave the risk set only after the
        // threshold's events (the events-before-censorings tie convention).
        while data_index < observations.len() && observations[data_index].y == threshold {
            let observation = observations[data_index];
            debug_assert!(!observation.observed);

            at_risk[observation.x] *=
                1.0 - observation.weight / covariate_statistics[observation.x].weight;
            covariate_statistics[observation.x].weight -= observation.weight;
            if data_index + 1 == last_positive[observation.x] {
                // Fully censored away: the remaining at-risk mass is exactly 0.
                at_risk[observation.x] = 0.0;
            }
            data_index += 1;
        }

        // Save results
        cdfs.extend(repeat_n(1.0, zero_count));
        let mut start_index = zero_count;
        for partition in partitions.drain(..) {
            // A pooled ratio target can exceed 1 (a group's raw Kaplan-Meier
            // survival sits above its pooled fitted survival after a neighbouring
            // group died, and zero-at-risk blocks don't pull it back). Survival
            // never increases, so cap the applied step at 1 — the exact solution
            // of the ratio regression under its θ ≤ 1 bound.
            let step = partition.value.min(1.0);
            for index in start_index..partition.index {
                // Previous iteration update
                at_risk[index] *= step;
                // Current iteration update
                survival[index] *= step;
                // Write out result, clamp to ensure numerical noise doesn't get us out of
                // [0, 1]
                cdfs.push(1.0 - survival[index].clamp(0.0, 1.0));
            }
            start_index = partition.index;
        }
        progress.increment();
    }

    cdfs
}

#[cfg(test)]
mod test {
    use crate::structures::Increasing;
    use crate::test::is_relative_eq_vec;
    use crate::total_order::hazard_rate_order::censored::algorithm;
    use crate::total_order::preprocessing::preprocess_censored;

    fn execute_test<const N: usize, const N_COVARIATE: usize, const N_THRESHOLD: usize>(
        x: [f64; N],
        y: [f64; N],
        observed: [bool; N],
        weight: [f64; N],
        expected: [[f64; N_COVARIATE]; N_THRESHOLD],
    ) {
        let context = preprocess_censored(&x, &y, &observed, &weight).unwrap();
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.unique_covariates.len(), N_COVARIATE);
        assert_eq!(context.thresholds.len(), N_THRESHOLD);
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

    /// An event at the covariate the PAVA cursor points to: each covariate is
    /// pushed exactly once per row. All observed, no pooling required, so the rows
    /// are plain conditional-survival bookkeeping.
    #[test]
    fn test_event_at_cursor() {
        execute_test(
            [0.0, 0.0, 0.0, 1.0],
            [1.0, 2.0, 9.0, 9.0],
            [true; 4],
            [1.0; 4],
            [[1.0 / 3.0, 0.0], [2.0 / 3.0, 0.0], [1.0, 1.0]],
        );
    }

    /// Censored-only response values are not thresholds: preprocessing snaps the
    /// censored observation onto the threshold below it, so the single-covariate
    /// Kaplan-Meier fast path emits exactly one value per threshold.
    /// KM: F(1) = 1/3; censor@2 leaves the at-risk set; F(3) = 1.
    #[test]
    fn test_censored_only_threshold_single_covariate() {
        execute_test(
            [5.0, 5.0, 5.0],
            [1.0, 2.0, 3.0],
            [true, false, true],
            [1.0; 3],
            [[1.0 / 3.0], [1.0]],
        );
    }

    /// Censored-only thresholds get no CDF row: the matrix is exactly
    /// `thresholds.len() x unique_covariates.len()`.
    #[test]
    fn test_censored_only_threshold_shape() {
        let context = preprocess_censored(
            &[0.0, 1.0, 2.0],
            &[1.0, 2.0, 3.0],
            &[true, false, true],
            &[1.0; 3],
        )
        .unwrap();
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(context.thresholds.len(), 2);
        assert_eq!(
            cdfs.len(),
            context.thresholds.len() * context.unique_covariates.len()
        );
        assert!(cdfs.iter().all(|v| (0.0..=1.0).contains(v)));
    }

    /// Hazard-rate order implies stochastic dominance: every row's CDF is
    /// nonincreasing in the covariate, including rows with several events where
    /// pooling must run after every partition push.
    #[test]
    fn test_covariate_ordering_invariant() {
        let context = preprocess_censored(
            &[2.0, 0.0, 1.0, 1.0, 1.0],
            &[2.0, 4.0, 1.0, 2.0, 1.0],
            &[true; 5],
            &[2.0, 2.0, 2.0, 1.0, 1.0],
        )
        .unwrap();
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        let n_cov = context.unique_covariates.len();
        assert_eq!(cdfs.len(), context.thresholds.len() * n_cov);
        for (t, row) in cdfs.chunks_exact(n_cov).enumerate() {
            for c in 1..n_cov {
                assert!(
                    row[c] <= row[c - 1] + 1e-6,
                    "threshold {t}: CDF increases with the covariate: {row:?}"
                );
            }
        }
    }

    /// Groups whose remaining mass is fully censored away (at-risk weight 0) still
    /// yield finite CDF values in [0, 1].
    #[test]
    fn test_zero_at_risk_groups() {
        let context = preprocess_censored(
            &[0.0, 1.0, 2.0, 2.0, 3.0, 3.0],
            &[9.0, 1.0, 2.0, 2.0, 1.0, 3.0],
            &[true, false, true, false, true, true],
            &[1.0; 6],
        )
        .unwrap();
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(
            cdfs.len(),
            context.thresholds.len() * context.unique_covariates.len()
        );
        assert!(
            cdfs.iter()
                .all(|v| v.is_finite() && (0.0..=1.0).contains(v))
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

    /// Fully observed data must match the uncensored kernel (same sequential
    /// least-squares problem, ratio-regression weights w·S(t-)²): the expected
    /// values are the uncensored `test_mixed_tonicity_3` / doctest values.
    #[test]
    fn test_mixed_tonicity_3() {
        execute_test(
            [0.0, 1.0, 2.0],
            [3.0, 1.0, 2.0],
            [true; 3],
            [1.0; 3],
            [
                [0.5, 0.5, 0.0],
                [5.0 / 6.0, 5.0 / 6.0, 2.0 / 3.0],
                [1.0, 1.0, 1.0],
            ],
        );
    }

    /// Fully observed data must match the uncensored kernel: the expected values
    /// are the uncensored `test_mixed_tonicity_6` values.
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
                [
                    25.0 / 34.0,
                    25.0 / 34.0,
                    11.0 / 17.0,
                    1.0 / 2.0,
                    1.0 / 2.0,
                    0.0,
                ],
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
                [1.0; 6],
            ],
        );
    }

    /// The covariates at 2.0 and 3.0 are censored strictly before the only event:
    /// they leave the risk set before every threshold, carry no information, and are
    /// discarded by preprocessing — the grid keeps [1.0, 4.0, 5.0].
    #[test]
    fn test_5() {
        execute_test(
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [3.0, 2.0, 1.0, 5.0, 4.0],
            [true, false, false, false, false],
            [1.0; 5],
            [[1.0, 0.0, 0.0]],
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

    /// A pooled ratio can exceed 1 in exact arithmetic, so the θ ≤ 1 cap on the
    /// applied step is part of the estimator, not noise cleanup. At t=1 the
    /// hazard 10/30 pools across both groups (fitted survival [2/3, 2/3] while
    /// the first group's raw Kaplan-Meier survival is still 1); the second
    /// group's remaining mass is then fully censored away (at-risk mass exactly
    /// 0, estimator frozen at 1/2). At t=2 the first group's raw ratio is
    /// 0.9 / (2/3) = 27/20 with positive weight, the frozen group's is 3/4 with
    /// weight 0 — pooling cannot pull the block back under 1. The cap clips the
    /// step to 1, so the t=2 row repeats the t=1 row; unclipped, survival would
    /// rise to 0.9 and the CDF would fall from 1/3 to 0.1, no longer a
    /// distribution. (The uncensored kernel never needs this: with constant
    /// weights the pooled ratios stay ≤ 1 in exact arithmetic; only censoring's
    /// non-uniform at-risk shrinkage breaks that bound.)
    #[test]
    fn test_frozen_group_ratio_above_one() {
        execute_test(
            [0.0, 0.0, 1.0, 1.0],
            [2.0, 3.0, 1.0, 1.0],
            [true, true, true, false],
            [1.0, 9.0, 10.0, 10.0],
            [[1.0 / 3.0, 1.0 / 3.0], [1.0 / 3.0, 1.0 / 3.0], [1.0, 0.5]],
        );
    }
}
