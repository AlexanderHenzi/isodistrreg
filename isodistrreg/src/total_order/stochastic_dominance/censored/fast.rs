use crate::progress::ProgressTracker;
use crate::routines::kaplan_meier;
use crate::structures::{Direction, Observation};
use crate::total_order;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use crate::total_order::stochastic_dominance::censored::propagate_bounds::{
    Avx2Kernel, Avx512Kernel,
};
use crate::total_order::stochastic_dominance::censored::propagate_bounds::{Kernel, ScalarKernel};
use crate::total_order::stochastic_dominance::censored::structures::{
    CompletionIndex, Estimates, Partition,
};
use crate::total_order::stochastic_dominance::routines;
use crate::total_order::structures;
use crate::total_order::structures::CensoredContext;
use crate::total_order::structures::WeightedPartition;
use crate::total_order::weight_noise_floor;
use std::cmp::Ordering;
use std::iter::repeat_n;

pub fn algorithm<D: Direction, X: crate::Float, Y: crate::Float>(
    context: &CensoredContext<X, Y>,
    progress: &dyn ProgressTracker,
) -> Vec<f32> {
    progress.set_total(context.n_threshold());

    if context.n_threshold() == 0 {
        debug_assert!(context.unique_covariates.is_empty());
        debug_assert!(context.thresholds.is_empty());
        return Vec::with_capacity(0);
    } else if context.n_threshold() == 1 {
        // Single threshold -> simple binary isotonic regression with censoring amount
        return total_order::routines::single_response::<D, _>(
            context.observations.clone(),
            context.covariate_statistics.clone(),
        );
    }
    if context.n_covariate() == 1 {
        // Single covariate -> a single empirical cdf
        return kaplan_meier(
            context.observations.iter().copied(),
            context.covariate_statistics[0].weight,
        );
    }

    // Collects final estimate
    let mut cdfs = Vec::with_capacity(context.n_threshold() * context.n_covariate());

    // Tracks which observation we're treating, sorted by response and censoring (and covariate)
    let mut data_index = 0;

    if data_index == context.n() {
        // No uncensored observations, we're done
        return cdfs;
    }

    // The next observation is uncensored
    debug_assert!(context.observations[data_index].observed);

    let just_one_more = data_index + 1 == context.n();
    if just_one_more {
        // Ends with a single uncensored observation, we're done
        finalize_for_single_uncensored::<D::REVERSE, _, _>(data_index, context, &mut cdfs);
        return cdfs;
    }

    // Build the exact-0 completion oracle once (O(n + C)); `Estimates` carries it
    // through the whole hot call tree so every cell write pins mathematically-completed
    // intervals to exactly 0 the moment they complete. None of the early returns above
    // need it: `single_response` already pins its completed covariates (`share = 1.0`)
    // and all-1.0 blocks never merge under the strict pooling comparison, `kaplan_meier`
    // tracks `last_positive_observed`, and `finalize_for_single_uncensored` emits
    // literal 0.0/1.0.
    let completion = CompletionIndex::new(&context.observations, context.n_covariate());

    let at_least_two_more = data_index + 1 < context.n();
    let at_least_two_uncensored = context.observations[data_index + 1].observed;
    let (start_threshold, estimates, partitions) = if at_least_two_more && at_least_two_uncensored {
        // Run the classical PAV algorithm first if we can save at least one of the more expensive
        // update steps of the censored algorithm.
        //
        // The uncensored prefix runs in f32 (classical PAVA — well-conditioned and shares helpers
        // with other algorithms). The bridge into the generalized-PAVA hot path happens at
        // `Estimates::from_partial_uncensored_solution`.

        // Buffer to temporarily store partitions right of the covariate index being changed
        let mut partitions_to_store = Vec::with_capacity(context.n_covariate());

        // Apply the algorithm for the uncensored case
        let (consumed_share, consumed_weight, partitions) = routines::accelerated_pava::<_, _, _, D>(
            &mut data_index,
            &context.observations,
            &context.covariate_statistics,
            |observed| !observed,
            &mut partitions_to_store,
            &mut cdfs,
            progress,
        );

        let final_threshold = context.observations.last().unwrap().y;
        if context.observations[data_index].y == final_threshold {
            // We're almost done, specialized final threshold and return
            if context.observations[data_index..]
                .iter()
                .all(|observation| !observation.observed)
            {
                // Fully uncensored final threshold
                cdfs.extend(repeat_n(1.0, context.n_covariate()));
            } else {
                // Some censoring in final threshold only
                finalize_for_censoring_only_in_final_threshold::<D, _, _>(
                    data_index,
                    consumed_share,
                    consumed_weight,
                    partitions,
                    context,
                    &mut cdfs,
                    partitions_to_store,
                    &completion,
                );
            }
            progress.increment();
            return cdfs;
        }

        // Initialize for the general algorithm
        let start_threshold = context.observations[data_index].y;
        let estimates = Estimates::from_partial_uncensored_solution(
            consumed_weight.weight,
            context.covariate_statistics.as_slice(),
            data_index,
            completion,
        );
        let index_only_partitions: Vec<_> = partitions
            .into_iter()
            .map(|structures::Partition { index, .. }| Partition::new(index))
            .collect();
        (start_threshold, estimates, index_only_partitions)
    } else {
        // Initialize directly for the more general algorithm

        // Set the start count already to the value that will be appropriate after initialization
        let mut estimates = Estimates::new(context.n_covariate(), data_index, completion);
        let mut partitions = Vec::with_capacity(context.n_covariate());

        // First uncensored threshold (fast initialization only)
        debug_assert!(context.observations[data_index].observed);
        initialize::<D::REVERSE, _, _>(data_index, &mut estimates, &mut partitions, context);
        let start_threshold = context.observations[data_index].y;
        data_index += 1;

        (start_threshold, estimates, partitions)
    };

    // Remaining thresholds with the full algorithm. The runtime CPU-feature check happens once
    // here; the chosen branch monomorphizes the entire inner call tree
    // (`generalized_pava → update_uncensored → pool → propagate_bounds*`) over the kernel's
    // zero-sized `Kernel` impl, so `K::apply` inlines into every callsite — no indirection in
    // the hot loop.
    dispatch_generalized_pava::<D, _, _>(
        data_index,
        start_threshold,
        estimates,
        partitions,
        context,
        &mut cdfs,
        progress,
    );

    cdfs
}

fn finalize_for_single_uncensored<D: Direction, X: crate::Float, Y: crate::Float>(
    data_index: usize,
    input: &CensoredContext<X, Y>,
    cdf: &mut Vec<f32>,
) {
    match D::IS_INCREASING {
        true => {
            let zeros_count = input.observations[data_index].x;
            cdf.extend(repeat_n(0.0, zeros_count));
            cdf.extend(repeat_n(1.0, input.n_covariate() - zeros_count));
        }
        false => {
            let ones_count = input.observations[data_index].x + 1;
            cdf.extend(repeat_n(1.0, ones_count));
            cdf.extend(repeat_n(0.0, input.n_covariate() - ones_count));
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn finalize_for_censoring_only_in_final_threshold<
    D: Direction,
    X: crate::Float,
    Y: crate::Float,
>(
    data_index: usize,
    mut consumed_share: Vec<f32>,
    mut consumed_weight: routines::ConsumedMass,
    mut partitions: Vec<WeightedPartition>,
    input: &CensoredContext<X, Y>,
    cdf: &mut Vec<f32>,
    mut partitions_to_store: Vec<WeightedPartition>,
    completion: &CompletionIndex,
) {
    for observation in &input.observations[data_index..] {
        if observation.observed {
            routines::classical_pava_update_step::<_, _, D::REVERSE>(
                observation,
                &mut consumed_share,
                &mut consumed_weight,
                &mut partitions,
                &input.covariate_statistics,
                &mut partitions_to_store,
            );
        }
    }

    // Pin completed partitions to exactly 1.0. Partition values here are CDF-space
    // shares (consumed/total of f32 running sums), which can sit ulps off 1.0 even when
    // every observation of the range was consumed as an event. In this path every
    // censored observation sorts at/after the final threshold's events, so with all data
    // consumed a range completes iff it contains no censoring at all — exactly the cells
    // whose CDF is mathematically 1. The partition→range mapping mirrors
    // `routines::store_in_cdf`: ascending indices for increasing fits (partition `l`
    // covers `[partitions[l - 1].index, partitions[l].index)`), descending otherwise.
    for l in 0..partitions.len() {
        let start = if D::IS_INCREASING {
            if l == 0 { 0 } else { partitions[l - 1].index }
        } else {
            partitions.get(l + 1).map_or(0, |p| p.index)
        };
        let end = partitions[l].index; // exclusive
        if completion.completes_with_all_data(start, end - 1) {
            partitions[l].value = 1.0;
        }
    }

    routines::store_in_cdf::<_, D>(&partitions, cdf);
}

fn initialize<D: Direction, X: crate::Float, Y: crate::Float>(
    data_index: usize,
    estimators: &mut Estimates,
    partitions: &mut Vec<Partition>,
    input: &CensoredContext<X, Y>,
) {
    let observation = &input.observations[data_index];
    let obs_x = observation.x;
    let obs_weight = observation.weight;

    // Initialize estimators. The loop runs r descending so the completion queries can
    // fold `marker(r)` into a running range-max over `[r, obs_x]`; the body is
    // order-independent (iteration r writes only cell `(r, obs_x)` plus the
    // `(r - 1, obs_x)` lower bound, and reads nothing written by other iterations).
    let mut marker_max = 0;
    for r in (0..=obs_x).rev() {
        marker_max = marker_max.max(estimators.completion.marker(r));
        debug_assert_eq!(marker_max, estimators.completion.range_max(r, obs_x));
        let total_weight = if r > 0 {
            input.covariate_statistics[obs_x].cumulative_weight
                - input.covariate_statistics[r - 1].cumulative_weight
        } else {
            input.covariate_statistics[obs_x].cumulative_weight
        };
        // A completing interval must be exactly 0 despite the f32 arithmetic: the cell's
        // total weight is a cumulative-weight difference, so even `1 - w/total` for a
        // lone observation is not exactly 0. Only `r == obs_x` can actually fire here
        // (any wider range contains another covariate whose observations all sort after
        // this very first one).
        let raw_value = if CompletionIndex::completes_marker(marker_max, data_index) {
            0.0
        } else {
            1.0 - obs_weight / total_weight
        };

        let (value, cold) = estimators.entry_mut(r, obs_x);
        cold.raw_value = raw_value;
        cold.weight = obs_weight;
        // The estimators (r, cov_index) are decreasing in r, so the lower bound is below the value
        *value = raw_value;
        cold.count = data_index + 1;

        // Propagate bound
        if r > 0 {
            estimators.cold_mut(r - 1, obs_x).lower_bound = raw_value;
        }
    }

    // Set up partitions
    match D::FORBIDDEN_ORDERING {
        Ordering::Less => {
            // The first antitonic regression has at least one partition
            partitions.push(Partition::new(obs_x + 1)); // Partition indices are exclusive
            // The first antitonic regression may have a second partition, if the value isn't the last
            if obs_x < input.n_covariate() - 1 {
                partitions.push(Partition::new(input.n_covariate()));
            }
        }
        Ordering::Greater => {
            if obs_x == 0 {
                partitions.push(Partition::new(input.n_covariate()));
            } else {
                partitions.push(Partition::new(obs_x));
                partitions.push(Partition::new(input.n_covariate()));
            }
        }
        Ordering::Equal => panic!(),
    }
}

/// Pick the kernel monomorphization based on runtime CPU features. AVX-512F → AVX2 → scalar.
fn dispatch_generalized_pava<D: Direction, X: crate::Float, Y: crate::Float>(
    data_index: usize,
    start_threshold: usize,
    estimates: Estimates,
    partitions: Vec<Partition>,
    input: &CensoredContext<X, Y>,
    cdf: &mut Vec<f32>,
    progress: &dyn ProgressTracker,
) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx512f") {
            return generalized_pava::<D, Avx512Kernel, _, _>(
                data_index,
                start_threshold,
                estimates,
                partitions,
                input,
                cdf,
                progress,
            );
        } else if is_x86_feature_detected!("avx2") {
            return generalized_pava::<D, Avx2Kernel, _, _>(
                data_index,
                start_threshold,
                estimates,
                partitions,
                input,
                cdf,
                progress,
            );
        }
    }

    generalized_pava::<D, ScalarKernel, _, _>(
        data_index,
        start_threshold,
        estimates,
        partitions,
        input,
        cdf,
        progress,
    );
}

fn generalized_pava<D: Direction, K: Kernel, X: crate::Float, Y: crate::Float>(
    mut data_index: usize,
    start_threshold: usize,
    mut estimates: Estimates,
    mut partitions: Vec<Partition>,
    input: &CensoredContext<X, Y>,
    cdf: &mut Vec<f32>,
    progress: &dyn ProgressTracker,
) {
    let mut tmp_partition_store = Vec::with_capacity(input.n_covariate());
    // Scratch buffers reused across all `pool` calls: `row_buf` holds `values[(r, k)]` for the
    // current outer-`r` iteration so the inner `s` sweep reads a contiguous slice instead of
    // striding into the triangle; `marker_buf` holds the completion-marker prefix maxima of the
    // current pooling round's ultimate range.
    let mut row_buf: Vec<f32> = Vec::with_capacity(input.n_covariate());
    let mut marker_buf: Vec<u32> = Vec::with_capacity(input.n_covariate());
    // Precompute the dynamic K-M safety threshold. In `update_value`, the K-M numerator
    // divides `obs.weight` by `remaining_weight = total_weight − cold.weight`. Both
    // `total_weight` and `cold.weight` are sums of f32 weights and pick up O(Σw · u_32)
    // absolute round-off, so `remaining_weight` can be non-zero by O(Σw · u_32) even when
    // every observation in the cell has already been applied. Without a guard the K-M
    // step would then divide by a sub-noise-floor value and blow up — see
    // `weight_noise_floor` for the bound's derivation.
    let epsilon = weight_noise_floor(input.observations.iter().map(|o| o.weight).sum());
    for threshold in start_threshold..input.n_threshold() {
        while data_index < input.n() {
            let observation = &input.observations[data_index];
            if observation.y > threshold {
                break;
            }

            if observation.observed {
                update_uncensored::<D, K, _, _>(
                    data_index,
                    &mut estimates,
                    &mut partitions,
                    input,
                    epsilon,
                    &mut tmp_partition_store,
                    &mut row_buf,
                    &mut marker_buf,
                );
            } else {
                // Censored observations are deferred. They affect the K-M estimate only at the
                // next uncensored arrival, which `update_value` picks up by walking forward
                // through `observations` from `self.count`.
            }

            data_index += 1;
        }

        store_in_cdf(&estimates, &partitions, cdf);
        progress.increment();
    }
}

#[allow(clippy::too_many_arguments)]
fn update_uncensored<D: Direction, K: Kernel, X: crate::Float, Y: crate::Float>(
    data_index: usize,
    estimates: &mut Estimates,
    partitions: &mut Vec<Partition>,
    input: &CensoredContext<X, Y>,
    epsilon: f32,
    tmp_partition_store: &mut Vec<Partition>,
    row_buf: &mut Vec<f32>,
    marker_buf: &mut Vec<u32>,
) {
    let observation = input.observations[data_index];
    let (partition_index, (lower, upper)) =
        routines::find_partition_bounds::<_, _, D::REVERSE>(observation.x, partitions);
    debug_assert!(lower <= observation.x && observation.x < upper);

    // Store right part of partitions
    tmp_partition_store.extend(partitions.drain(partition_index + 1..));

    if D::FORBIDDEN_ORDERING == Ordering::Less {
        // TODO
        unimplemented!("Need to reverse the partition management and pooling");
    }

    // Split the triangle and update computation of left sub-triangle
    partitions[partition_index].index = observation.x + 1; // partition indices are exclusive
    // Update the triangle of the range of the new partition
    estimates.update_partial_row_with_single_observation::<K, _, _>(
        data_index,
        lower,
        &observation,
        input,
        epsilon,
    );
    // Pooling left part of partitions (direction is the same, because we're working with survival
    // quantities, not the CDF)
    pool::<_, _, D, K, _, _>(
        data_index,
        observation.x,
        estimates,
        partitions,
        input,
        epsilon,
        row_buf,
        marker_buf,
    );

    // Accelerated extension and pooling
    for i in observation.x + 1..upper {
        partitions.push(Partition::new(i + 1));
        // Direction is the same, because we're working with survival quantities, not the CDF)
        pool::<_, _, D, K, _, _>(
            data_index,
            observation.x,
            estimates,
            partitions,
            input,
            epsilon,
            row_buf,
            marker_buf,
        );
    }

    // Restore right-most partitions
    partitions.append(tmp_partition_store);
}

fn pool<W, V, D: Direction, K: Kernel, X: crate::Float, Y: crate::Float>(
    data_index: usize,
    split_x: usize,
    estimates: &mut Estimates,
    partitions: &mut Vec<structures::Partition<W, V>>,
    input: &CensoredContext<X, Y>,
    epsilon: f32,
    row_buf: &mut Vec<f32>,
    marker_buf: &mut Vec<u32>,
) {
    loop {
        // Start inclusive, end exclusive
        let (penultimate_start, penultimate_end, ultimate_start, ultimate_end) = match partitions
            .as_slice()
        {
            [.., antepenultimate, penultimate, ultimate] => (
                antepenultimate.index,
                penultimate.index,
                penultimate.index,
                ultimate.index,
            ),
            [penultimate, ultimate] => (0, penultimate.index, penultimate.index, ultimate.index),
            _ => break,
        };

        let penultimate_value = estimates.value(penultimate_start, penultimate_end - 1);
        let ultimate_value = estimates.value(ultimate_start, ultimate_end - 1);

        if penultimate_value.partial_cmp(&ultimate_value).unwrap() != D::FORBIDDEN_ORDERING {
            break;
        }
        // Perform pooling
        let ultimate_boundary = partitions.pop().unwrap();
        *partitions.last_mut().unwrap() = ultimate_boundary;

        if penultimate_start > split_x {
            // Both merged blocks lie strictly right of the covariate that triggered this
            // update, so every cross cell (r, s) satisfies split_x < r <= s. The new
            // observation at split_x is outside [r, s], and both blocks are sub-intervals
            // of the block that was just split — whose interior cells were all current
            // before the split (censored arrivals never change stored values, and an
            // uncensored arrival inside would have split it earlier). Values, bounds, and
            // clips are therefore already exactly what the loop below would recompute,
            // and the completion state cannot newly fire either: a completing event with
            // x in [r, s] and index <= data_index would itself have split the block.
            // Keep the merge, skip the no-op cell updates.
            continue;
        }

        // The completion range-max over a cell `[r, s]` splits at the partition boundary
        // (`ultimate_start == penultimate_end` in both match arms above): the left part
        // `[r, penultimate_end)` extends one covariate per descending `r` and folds
        // incrementally below; the right part `[ultimate_start, s]` is shared by every
        // `r`, so its prefix maxima are materialized once per pooling round —
        // `marker_buf[s - ultimate_start]` = max of `m[ultimate_start..=s]`. This keeps
        // the per-cell query O(1) without touching the row-materialization loop.
        marker_buf.clear();
        let mut ultimate_marker_max = 0;
        for x in ultimate_start..ultimate_end {
            ultimate_marker_max = ultimate_marker_max.max(estimates.completion.marker(x));
            marker_buf.push(ultimate_marker_max);
        }
        let mut left_marker_max = 0;

        // TODO: Can we iterate in a different order that tries to find a min == max combination
        //  as soon as possible?
        // Iterating over r in the outer loop and over s in the inner loop is somehow significantly
        // faster
        for r in (penultimate_start..penultimate_end).rev() {
            left_marker_max = left_marker_max.max(estimates.completion.marker(r));
            // Materialize row r once: `row_buf[k - r] = values[(r, k)]` for `k` in
            // `r..right_partition_end`. The inner `s` sweep reuses this row across every `s`, and
            // refreshes the just-written entry after `update_value`.
            row_buf.clear();
            row_buf.reserve(ultimate_end - r);
            let mut col_k_base = r * (r + 1) / 2;
            row_buf.push(estimates.values[col_k_base + r]); // (r, r)
            for k in (r + 1)..penultimate_end {
                col_k_base += k;
                row_buf.push(estimates.values[col_k_base + r]);
            }

            for s in ultimate_start..ultimate_end {
                let marker_max = left_marker_max.max(marker_buf[s - ultimate_start]);
                estimates.propagate_bounds_with_row::<K>(r, s, row_buf);
                estimates.update_value(data_index, r, s, input, epsilon, marker_max);
                // Reflect the freshly written `values[(r, s)]` back into the row buffer so
                // subsequent iterations (with larger `s`) see the updated value.
                let idx = Estimates::compute_index((r, s), estimates.len());
                row_buf.push(estimates.values[idx]);
            }
        }
    }
}

impl Estimates {
    fn update_value<X: crate::Float, Y: crate::Float>(
        &mut self,
        data_index: usize,
        covariate_start_index: usize, // inclusive
        covariate_end_index: usize,   // inclusive
        input: &CensoredContext<X, Y>,
        epsilon: f32,
        // Caller-maintained running max of `completion.marker` over the cell's covariate
        // range — both callers extend their interval one covariate at a time, which makes
        // the range-max O(1) without any index structure.
        marker_max: u32,
    ) {
        debug_assert_eq!(
            marker_max,
            self.completion
                .range_max(covariate_start_index, covariate_end_index),
        );
        // Exact-0 pinning, checked before anything else: once the interval's last
        // observation has been consumed (an event — see `CompletionIndex`), the interval
        // Kaplan-Meier survival is mathematically exactly 0 and must be stored as exactly
        // 0 despite f32 drift. Deliberately UNCLIPPED: bounds freshly propagated from
        // not-yet-exact neighbors can sit ulps above 0 and must not lift the exact zero.
        // Deliberately no `cold.count`/`cold.weight` update: no in-range observation
        // remains unconsumed, and this branch re-fires on every future call (the
        // condition is monotone in `data_index`), so the stale bookkeeping is never read
        // again.
        if CompletionIndex::completes_marker(marker_max, data_index) {
            let (value, cold) = self.entry_mut(covariate_start_index, covariate_end_index);
            cold.raw_value = 0.0;
            *value = 0.0;
            return;
        }

        let (value, cold) = self.entry_mut(covariate_start_index, covariate_end_index);

        // `cold.count` may already be `data_index + 1` when the cell was touched earlier in the
        // same run of tied uncensored observations (the walk consumes the whole run at once).
        debug_assert!(data_index + 1 >= cold.count);
        debug_assert!(covariate_start_index <= covariate_end_index);
        debug_assert!(covariate_end_index < input.n_covariate());

        // We assume the bounds to be up to date. If a bound is not available, will be false.
        debug_assert_eq!(
            cold.lower_bound.is_nan(),
            cold.upper_bound.is_nan(),
            "Bounds should only be NAN if this is a diagonal value",
        );
        if cold.lower_bound >= cold.upper_bound {
            *value = f32::midpoint(cold.lower_bound, cold.upper_bound);
            return;
        }

        let total_weight = if covariate_start_index > 0 {
            input.covariate_statistics[covariate_end_index].cumulative_weight
                - input.covariate_statistics[covariate_start_index - 1].cumulative_weight
        } else {
            input.covariate_statistics[covariate_end_index].cumulative_weight
        };
        let mut remaining_weight = total_weight - cold.weight;
        if remaining_weight < epsilon {
            // There is no more weight to process, so the raw Kaplan-Meier value is final — but
            // the bounds may have just been re-propagated from fresh neighbor values, so the
            // clip still has to be reapplied. Skipping it would leave a value clipped against
            // bounds that no longer exist.
            *value = cold.raw_value.max(cold.lower_bound).min(cold.upper_bound);
            return;
        }

        let mut raw_value = cold.raw_value;

        // Walk the response-sorted `observations` slice forward from `cold.count`, up to and
        // including the observation currently being processed — and no further. Observations
        // are sorted by (response asc, uncensored first, covariate asc), so everything in this
        // range is at or below the current threshold and must be folded in. Reading *ahead* of
        // `data_index` (e.g. absorbing the rest of a run of tied uncensored responses) is not
        // allowed: `cold.count` could then no longer describe what was consumed, and the cell's
        // raw value would run ahead of the neighbor cells its bounds are computed from.
        // No sort or temp buffer is needed: items are already in K-M apply order.
        for obs in &input.observations[cold.count..=data_index] {
            if obs.x < covariate_start_index || obs.x > covariate_end_index {
                continue;
            }
            if obs.observed {
                raw_value *= 1.0 - obs.weight / remaining_weight;
            }
            remaining_weight -= obs.weight;
        }

        cold.raw_value = raw_value;
        cold.weight = total_weight - remaining_weight;
        cold.count = data_index + 1;

        *value = raw_value.max(cold.lower_bound).min(cold.upper_bound);
    }
}

impl Estimates {
    /// Update a row in the triangle.
    ///
    /// We need to update only a row (and not the entire triangle) because the items in the triangle
    /// not part of this row are assumed to be up to date because they were part of the previous
    /// partition and not affected by the new single observation.
    fn update_partial_row_with_single_observation<K: Kernel, X: crate::Float, Y: crate::Float>(
        &mut self,
        data_index: usize,
        partition_start_index: usize,
        observation: &Observation<usize, usize, bool, f32>,
        input: &CensoredContext<X, Y>,
        epsilon: f32,
    ) {
        // Completion range-max over `[r, observation.x]`, folded as the loop descends.
        let mut marker_max = 0;
        for r in (partition_start_index..=observation.x).rev() {
            marker_max = marker_max.max(self.completion.marker(r));
            // TODO: Try eliminating this branch
            if r < observation.x {
                self.propagate_bounds::<K>(r, observation.x);
            }
            self.update_value(data_index, r, observation.x, input, epsilon, marker_max);
        }
    }

    /// Propagate bounds — reads row r directly from the strided triangle. Used by the one-shot
    /// non-pool callsite (`update_partial_row_with_single_observation`) where there is no row
    /// reuse. Uses a stack array for small rows; heap fallback for large ones.
    fn propagate_bounds<K: Kernel>(&mut self, r: usize, s: usize) {
        assert!(r < s);
        assert!(s < self.len());

        let col_s_base = s * (s + 1) / 2;
        let col = &self.values[col_s_base + r + 1..=col_s_base + s];

        // Gather the row r entries (r, r), (r, r+1), ..., (r, s-1) into a small stack-bounded
        // local. The triangle is strided in the row direction, so we read scalars one at a
        // time; the kernel then sees two contiguous slices and can auto-vectorize.
        let len = s - r;
        let mut row_buf = [0.0f32; 128];
        let row_slice = if len <= 128 {
            let mut col_i_base = r * (r + 1) / 2;
            row_buf[0] = self.values[col_i_base + r];
            #[allow(clippy::needless_range_loop)]
            for k in 1..len {
                col_i_base += r + k;
                row_buf[k] = self.values[col_i_base + r];
            }
            &row_buf[..len]
        } else {
            // Fallback for unusually large rows: heap allocation instead of stack overflow.
            let mut v = Vec::with_capacity(len);
            let mut col_i_base = r * (r + 1) / 2;
            v.push(self.values[col_i_base + r]);
            for k in 1..len {
                col_i_base += r + k;
                v.push(self.values[col_i_base + r]);
            }
            let (lower, upper) = K::apply(&v, col);
            let cold = &mut self.cold[col_s_base + r];
            cold.lower_bound = lower;
            cold.upper_bound = upper;
            return;
        };

        let (lower, upper) = K::apply(row_slice, col);
        let cold = &mut self.cold[col_s_base + r];
        cold.lower_bound = lower;
        cold.upper_bound = upper;
    }

    /// Propagate bounds at `(r, s)` using a precomputed row buffer.
    ///
    /// `row_buf[k - r]` must equal `values[(r, k)]` for `k = r..s` (entries beyond `s` are
    /// ignored). Used in `pool` where the same row r is reused across an entire inner s sweep.
    fn propagate_bounds_with_row<K: Kernel>(&mut self, r: usize, s: usize, row_buf: &[f32]) {
        assert!(r < s);
        assert!(s < self.len());
        debug_assert!(row_buf.len() >= s - r);

        let len = s - r;
        let col_s_base = s * (s + 1) / 2;
        let row = &row_buf[..len];
        let col = &self.values[col_s_base + r + 1..=col_s_base + s];

        let (lower, upper) = K::apply(row, col);

        let cold = &mut self.cold[col_s_base + r];
        cold.lower_bound = lower;
        cold.upper_bound = upper;
    }
}

/// Store the latest antitonic regression as represented by the partitions in the CDF.
///
/// We clamp before saving, because we sometimes shave off a tiny epsilon above 1.0, for example.
fn store_in_cdf<W, V>(
    estimates: &Estimates,
    partitions: &[structures::Partition<W, V>],
    cdf: &mut Vec<f32>,
) {
    let partition_len = partitions[0].index;
    let value = estimates.value(0, partitions[0].index - 1);
    cdf.extend(repeat_n(1.0 - value.clamp(0.0, 1.0), partition_len));
    for l in 1..partitions.len() {
        let partition_len = partitions[l].index - partitions[l - 1].index;
        let value = estimates.value(partitions[l - 1].index, partitions[l].index - 1);
        cdf.extend(repeat_n(1.0 - value.clamp(0.0, 1.0), partition_len));
    }
}

#[cfg(test)]
mod test {
    use crate::structures::Increasing;
    use crate::total_order::preprocessing::preprocess_censored as preprocess;
    use crate::total_order::stochastic_dominance::censored::algorithm;

    #[test]
    fn small() {
        let context = preprocess(
            &[-29., -19., -18., -33., -23.],
            &[3., 23., 165., 5., 57.],
            &[false, false, true, false, true],
            &[1.0; 5],
        )
        .unwrap();
        let cdfs = algorithm::<Increasing, _, _>(&context, &crate::NoProgress);
        assert_eq!(cdfs, vec![1., 0., 1., 1.]);
        assert_eq!(context.unique_covariates, vec![-23., -18.]);
        assert_eq!(context.thresholds, vec![57., 165.]);
    }
}
