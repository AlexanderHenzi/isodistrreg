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
    Bounds, CompletionIndex, Estimates, Partition,
};
use crate::total_order::stochastic_dominance::routines;
use crate::total_order::structures;
use crate::total_order::structures::CensoredContext;
use crate::total_order::structures::WeightedPartition;
use crate::total_order::weight_noise_floor;
use std::cmp::Ordering;
use std::iter::repeat_n;

mod sweep;
use sweep::{SweepArgs, SweepDispatch};

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
            &context.observations,
        );
        let index_only_partitions: Vec<_> = partitions
            .into_iter()
            .map(|structures::Partition { index, .. }| Partition::new(index))
            .collect();
        (start_threshold, estimates, index_only_partitions)
    } else {
        // Initialize directly for the more general algorithm

        // Set the start count already to the value that will be appropriate after initialization
        let mut estimates = Estimates::new(
            context.n_covariate(),
            data_index,
            completion,
            &context.observations,
            &context.covariate_statistics,
        );
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
    // (`generalized_pava → update_threshold → pool → propagate_bounds*`) over the kernel's
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
            routines::update_threshold::<_, _, D::REVERSE>(
                std::slice::from_ref(observation),
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

        let idx = Estimates::compute_index((r, obs_x), estimators.len());
        let row_idx = Estimates::compute_row_index((r, obs_x), estimators.len());
        let cold = &mut estimators.cold[idx];
        cold.raw_value = raw_value;
        cold.weight = obs_weight;
        cold.count = (data_index + 1) as u32;
        estimators.set_value(idx, row_idx, raw_value);
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

fn generalized_pava<D: Direction, K: SweepDispatch, X: crate::Float, Y: crate::Float>(
    mut data_index: usize,
    start_threshold: usize,
    mut estimates: Estimates,
    mut partitions: Vec<Partition>,
    input: &CensoredContext<X, Y>,
    cdf: &mut Vec<f32>,
    progress: &dyn ProgressTracker,
) {
    let mut tmp_partition_store = Vec::with_capacity(input.n_covariate());
    // Scratch buffers reused across all `pool` calls: the completion-marker prefix maxima
    // of the current pooling round's ultimate range, and the suffix maxima of its
    // penultimate range. (The kernel's row operands come straight from the row-major
    // `values_row` mirror — no row scratch is needed.)
    let mut marker_buf: Vec<u32> = Vec::with_capacity(input.n_covariate());
    let mut left_marker_buf: Vec<u32> = Vec::with_capacity(input.n_covariate());
    // Precompute the dynamic K-M safety threshold. In `update_value`, the K-M numerator
    // divides `obs.weight` by `remaining_weight = total_weight − cold.weight`. Both
    // `total_weight` and `cold.weight` are sums of f32 weights and pick up O(Σw · u_32)
    // absolute round-off, so `remaining_weight` can be non-zero by O(Σw · u_32) even when
    // every observation in the cell has already been applied. Without a guard the K-M
    // step would then divide by a sub-noise-floor value and blow up — see
    // `weight_noise_floor` for the bound's derivation. The total is taken from the
    // cumulative covariate statistics — the same accumulation the guarded
    // `total_weight` values are differences of.
    let epsilon = weight_noise_floor(input.covariate_statistics.last().unwrap().cumulative_weight);
    // Each iteration treats exactly one threshold: its events form one contiguous run
    // (observations sort by response, events before censorings, covariate ascending;
    // thresholds are exactly the unique event responses, so every threshold owns such
    // a run), followed by the threshold's censored observations. Only the very first
    // iteration can instead begin mid-threshold: after the uncensored-prefix bridge
    // `data_index` sits at the censored remainder whose events the prefix already
    // consumed, and after `initialize` it sits at whatever follows the first event.
    for threshold in start_threshold..input.n_threshold() {
        if data_index < input.n() {
            let observation = &input.observations[data_index];
            if observation.observed && observation.y == threshold {
                // The estimator is defined per threshold, so the whole event run is
                // applied as one batched update: all folds run at the run's last data
                // index and every pooled pair is swept exactly once, where per-arrival
                // processing re-splits and re-sweeps overlapping rectangles once per
                // arrival (measured 2.6-4.3x cell-visit multiplicity on tied-response
                // inputs).
                let mut run_end = data_index + 1;
                while run_end < input.n() {
                    let next = &input.observations[run_end];
                    if !next.observed || next.y != threshold {
                        break;
                    }
                    run_end += 1;
                }
                update_threshold::<D, K, _, _>(
                    data_index,
                    run_end,
                    &mut estimates,
                    &mut partitions,
                    input,
                    epsilon,
                    &mut tmp_partition_store,
                    &mut marker_buf,
                    &mut left_marker_buf,
                );
                data_index = run_end;
            }
        }

        // The threshold's censored observations are deferred. They affect the K-M
        // estimate only at the next event, whose walk picks them up by walking
        // forward through `observations` from `self.count`.
        while data_index < input.n() && input.observations[data_index].y == threshold {
            debug_assert!(
                !input.observations[data_index].observed,
                "events precede censorings within a threshold",
            );
            data_index += 1;
        }

        store_in_cdf(&estimates, &partitions, cdf);
        progress.increment();
    }
}

/// Apply one threshold's events — the run of uncensored observations tied on that
/// response, `observations[run_start..run_end]`, covariates strictly ascending by the
/// sort order — as a single batched update. Censored observations never trigger
/// updates and a threshold's events are contiguous, so each call handles exactly one
/// threshold. A run of one is exactly the classical per-arrival update step —
/// split-inheritance, leftward pooling, and the trailing flush as the accelerated
/// extension — so there is no separate single-arrival path.
///
/// The estimator is defined per threshold (see `definition`), so only the state after
/// the whole run is ever read: `store_in_cdf` runs once per threshold, and every
/// pooling comparison made *during* the run may already use threshold-end values.
/// Batching therefore folds every cell to the run's last data index the first time it
/// is touched, and arranges the partition edits so that nothing built during the run is
/// ever split again within the run:
///
/// - Arrivals are processed in ascending covariate order. A block formed by pooling at
///   arrival `x_i` spans covariates `<= x_i < x_{i+1}`, so no later arrival of the run
///   splits it — every pooled pair is swept exactly once. Per-arrival processing
///   instead re-splits the freshly rebuilt right-hand structure at each next arrival
///   and re-sweeps the overlapping rectangles, once per arrival.
/// - The accelerated extension (re-adding a split block's right remainder as singleton
///   blocks) is deferred until the remainder is actually reached — by the next arrival
///   needing to pass it, or by the final restore. The pool calls for these re-adds pass
///   `split_x = x_prev` (the arrival that split the block): every cell with both
///   coordinates in the remainder lies in a sub-interval of that just-split block and
///   contains no processed arrival, so the whole-rectangle skip and row clamp apply
///   exactly as in the per-arrival flow.
/// - Untouched whole blocks between arrivals are restored without pooling, like the
///   caller's blind restore of `tmp_partition_store`: their values and their left
///   boundary's validity are unchanged. Blocks the leftward cascades must merge them
///   with are reached through arrival pools, whose rectangles sweep full penultimate
///   rows (their `split_x` is the arrival covariate, to their right).
/// - An arrival inside an untouched original block splits it by inheritance
///   (`[block_start, x + 1)` stays one block), the same move as the per-arrival flow.
///   An arrival inside a deferred remainder enters as a singleton block and the pool
///   comparisons decide all merging — the canonical PAVA insertion that the extension
///   already uses for every other remainder covariate.
///
/// In exact arithmetic the result equals per-arrival processing (both evaluate the
/// per-threshold definition); in f32 a pooling comparison between near-equal heads can
/// resolve differently because batched comparisons see threshold-end values where
/// per-arrival comparisons see mid-run folds. The differential suites against
/// `definition` gate both schedules the same way.
#[allow(clippy::too_many_arguments)]
fn update_threshold<D: Direction, K: SweepDispatch, X: crate::Float, Y: crate::Float>(
    run_start: usize,
    run_end: usize,
    estimates: &mut Estimates,
    partitions: &mut Vec<Partition>,
    input: &CensoredContext<X, Y>,
    epsilon: f32,
    tmp_partition_store: &mut Vec<Partition>,
    marker_buf: &mut Vec<u32>,
    left_marker_buf: &mut Vec<u32>,
) {
    // Threshold-end folds throughout: every walk consumes the whole run.
    let data_index = run_end - 1;

    if D::FORBIDDEN_ORDERING == Ordering::Less {
        // TODO
        unimplemented!("Need to reverse the partition management and pooling");
    }

    // First arrival: the per-arrival split-inheritance move on its containing block.
    let first = input.observations[run_start];
    let (partition_index, (lower, upper)) =
        routines::find_partition_bounds::<_, _, D::REVERSE>(first.x, partitions);
    debug_assert!(lower <= first.x && first.x < upper);
    tmp_partition_store.extend(partitions.drain(partition_index + 1..));
    partitions[partition_index].index = first.x + 1;
    estimates.update_partial_row_with_single_observation::<K, _, _>(
        data_index, lower, &first, input, epsilon,
    );
    pool::<_, _, D, K, _, _>(
        data_index,
        first.x,
        estimates,
        partitions,
        input,
        epsilon,
        marker_buf,
        left_marker_buf,
    );

    // Covariates `[cursor, pending_end)` are not on the stack yet. With
    // `pending_is_remainder` they are the deferred right remainder of the block that
    // arrival `x_prev` split (re-added as singletons when reached); otherwise they are
    // an untouched original block (restored whole). `tmp_partition_store[store_cursor..]`
    // holds the original blocks not yet reached, in covariate order.
    let mut cursor = first.x + 1;
    let mut pending_end = upper;
    let mut pending_is_remainder = true;
    let mut x_prev = first.x;
    let mut store_cursor = 0usize;

    for arrival in run_start + 1..run_end {
        let observation = input.observations[arrival];
        let x = observation.x;
        debug_assert!(x > x_prev, "tied uncensored covariates are deduplicated");

        // Advance the stack up to the pending region containing x.
        while pending_end <= x {
            if pending_is_remainder {
                for i in cursor..pending_end {
                    partitions.push(Partition::new(i + 1));
                    pool::<_, _, D, K, _, _>(
                        data_index,
                        x_prev,
                        estimates,
                        partitions,
                        input,
                        epsilon,
                        marker_buf,
                        left_marker_buf,
                    );
                }
            } else {
                partitions.push(Partition::new(pending_end));
            }
            cursor = pending_end;
            pending_end = tmp_partition_store[store_cursor].index;
            pending_is_remainder = false;
            store_cursor += 1;
        }

        if pending_is_remainder {
            // x lies in a deferred remainder: re-add the covariates left of it, then
            // insert x as its own singleton block; pooling decides all merging.
            for i in cursor..x {
                partitions.push(Partition::new(i + 1));
                pool::<_, _, D, K, _, _>(
                    data_index,
                    x_prev,
                    estimates,
                    partitions,
                    input,
                    epsilon,
                    marker_buf,
                    left_marker_buf,
                );
            }
            partitions.push(Partition::new(x + 1));
            estimates.update_partial_row_with_single_observation::<K, _, _>(
                data_index,
                x,
                &observation,
                input,
                epsilon,
            );
        } else {
            // x lies in an untouched original block `[cursor, pending_end)`:
            // split-inheritance, exactly like the first arrival.
            partitions.push(Partition::new(x + 1));
            estimates.update_partial_row_with_single_observation::<K, _, _>(
                data_index,
                cursor,
                &observation,
                input,
                epsilon,
            );
        }
        pool::<_, _, D, K, _, _>(
            data_index,
            x,
            estimates,
            partitions,
            input,
            epsilon,
            marker_buf,
            left_marker_buf,
        );
        cursor = x + 1;
        pending_is_remainder = true;
        x_prev = x;
    }

    // Trailing remainder of the last split block, then the blind restore of the
    // untouched right part — the same closing moves as the per-arrival flow.
    for i in cursor..pending_end {
        partitions.push(Partition::new(i + 1));
        pool::<_, _, D, K, _, _>(
            data_index,
            x_prev,
            estimates,
            partitions,
            input,
            epsilon,
            marker_buf,
            left_marker_buf,
        );
    }
    partitions.extend(tmp_partition_store.drain(store_cursor..));
    tmp_partition_store.clear();
}

#[allow(clippy::too_many_arguments)]
fn pool<W, V, D: Direction, K: SweepDispatch, X: crate::Float, Y: crate::Float>(
    data_index: usize,
    split_x: usize,
    estimates: &mut Estimates,
    partitions: &mut Vec<structures::Partition<W, V>>,
    input: &CensoredContext<X, Y>,
    epsilon: f32,
    marker_buf: &mut Vec<u32>,
    left_marker_buf: &mut Vec<u32>,
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
        // Suffix maxima of the penultimate range: `left_marker_buf[r - penultimate_start]`
        // = max of `m[r..penultimate_end]`, so any cell-visit schedule can read its left
        // part in O(1).
        left_marker_buf.clear();
        left_marker_buf.resize(penultimate_end - penultimate_start, 0);
        let mut left_marker_max = 0;
        for r in (penultimate_start..penultimate_end).rev() {
            left_marker_max = left_marker_max.max(estimates.completion.marker(r));
            left_marker_buf[r - penultimate_start] = left_marker_max;
        }

        // The block-level skip above applies row-wise too: a cell (r, s) with r > split_x
        // has split_x outside [r, s] and, because every partition in play here lies within
        // the just-split block's original range (the partitions right of it were drained
        // to `tmp_partition_store`), [r, s] is a sub-interval of that block — so by the
        // same argument its values, bounds, clips, and completion state are already
        // exactly what the sweep would recompute. Only rows r <= split_x can change; clamp
        // the sweep to them. (`row_end > penultimate_start` holds here — otherwise the
        // whole-rectangle skip above fired.)
        let row_end = penultimate_end.min(split_x + 1);
        debug_assert!(row_end > penultimate_start);

        K::sweep_rect(SweepArgs {
            estimates,
            data_index,
            penultimate_start,
            row_end,
            ultimate_start,
            ultimate_end,
            input,
            epsilon,
            marker_buf,
            left_marker_buf,
        });
    }
}

impl Estimates {
    /// Update the K-M value of cell `idx` = (covariate_start_index, covariate_end_index),
    /// clipping against `bounds` (the cell's current clip bounds, usually just computed by
    /// `propagate_bounds_with_row` and handed over in registers).
    #[allow(clippy::too_many_arguments)]
    fn update_value<K: Kernel, X: crate::Float, Y: crate::Float>(
        &mut self,
        data_index: usize,
        idx: usize,
        bounds: Bounds,
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
        debug_assert_eq!(
            idx,
            Self::compute_index((covariate_start_index, covariate_end_index), self.len()),
        );
        let row_idx =
            Self::compute_row_index((covariate_start_index, covariate_end_index), self.len());
        let (lower_bound, upper_bound) = bounds;
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
            debug_assert!(idx < self.cold.len());
            // SAFETY: `idx` is the triangle index of the cell (asserted equal to
            // `compute_index` above), hence within `cold`.
            unsafe { self.cold.get_unchecked_mut(idx) }.raw_value = 0.0;
            self.set_value(idx, row_idx, 0.0);
            return;
        }

        debug_assert!(covariate_start_index <= covariate_end_index);
        debug_assert!(covariate_end_index < input.n_covariate());

        // We assume the bounds to be up to date. If a bound is not available, will be false.
        debug_assert_eq!(
            lower_bound.is_nan(),
            upper_bound.is_nan(),
            "Bounds should only be NAN if this is a diagonal value",
        );
        if lower_bound >= upper_bound {
            self.set_value(idx, row_idx, f32::midpoint(lower_bound, upper_bound));
            return;
        }

        // Branch-free total weight off the flat shifted `cum_w` copy — bitwise equal to
        // the `covariate_start_index > 0` branching form (`cum_w[0]` is exactly 0.0 and
        // `a - 0.0 == a` in IEEE arithmetic).
        debug_assert_eq!(
            self.cum_w[covariate_end_index + 1] - self.cum_w[covariate_start_index],
            if covariate_start_index > 0 {
                input.covariate_statistics[covariate_end_index].cumulative_weight
                    - input.covariate_statistics[covariate_start_index - 1].cumulative_weight
            } else {
                input.covariate_statistics[covariate_end_index].cumulative_weight
            }
        );
        // SAFETY: `covariate_end_index < self.len()` (asserted above) and `cum_w` has
        // `self.len() + 1` entries, so both indices are in bounds.
        let total_weight = unsafe {
            self.cum_w.get_unchecked(covariate_end_index + 1)
                - self.cum_w.get_unchecked(covariate_start_index)
        };

        debug_assert!(idx < self.cold.len());
        // SAFETY: `idx` is the triangle index of the cell (asserted equal to
        // `compute_index` above), hence within `cold`.
        let cold = unsafe { self.cold.get_unchecked_mut(idx) };
        // `cold.count` may already be `data_index + 1` when the cell was touched earlier in the
        // same run of tied uncensored observations (the walk consumes the whole run at once).
        debug_assert!(data_index + 1 >= cold.count as usize);
        let remaining_weight = total_weight - cold.weight;
        if remaining_weight < epsilon {
            // There is no more weight to process, so the raw Kaplan-Meier value is final — but
            // the bounds may have just been re-propagated from fresh neighbor values, so the
            // clip still has to be reapplied. Skipping it would leave a value clipped against
            // bounds that no longer exist.
            let result = cold.raw_value.max(lower_bound).min(upper_bound);
            self.set_value(idx, row_idx, result);
            return;
        }

        // Walk the response-sorted packed walk arrays forward from `cold.count`, up to and
        // including the observation currently being processed — and no further. Observations
        // are sorted by (response asc, uncensored first, covariate asc), so everything in this
        // range is at or below the current threshold and must be folded in. Reading *ahead* of
        // `data_index` (e.g. absorbing the rest of a run of tied uncensored responses) is not
        // allowed: `cold.count` could then no longer describe what was consumed, and the cell's
        // raw value would run ahead of the neighbor cells its bounds are computed from.
        // No sort or temp buffer is needed: items are already in K-M apply order.
        //
        // The filter reads only the 4-byte `walk_x` stream; the weight is loaded just for
        // items that pass. `K::walk_scan` vectorizes the range test where SIMD masks are
        // available — the walk is usually filter-dominated, with only a few percent of
        // items passing.
        let from = cold.count as usize;
        let to = data_index + 1;
        debug_assert!(from <= to && to <= self.walk_x.len() && to <= self.walk_w.len());
        // SAFETY: `from <= to` (asserted above from `cold.count`'s invariant) and
        // `to == data_index + 1 <= n`, the length of both walk arrays.
        let (walk_xs, walk_ws) = unsafe {
            (
                self.walk_x.get_unchecked(from..to),
                self.walk_w.get_unchecked(from..to),
            )
        };
        let (raw_value, remaining_weight) = K::walk_scan(
            walk_xs,
            walk_ws,
            covariate_start_index as u32,
            (covariate_end_index - covariate_start_index) as u32,
            cold.raw_value,
            remaining_weight,
        );

        // SAFETY: as above; re-borrowed because the walk's slice borrows ended the
        // earlier exclusive borrow.
        let cold = unsafe { self.cold.get_unchecked_mut(idx) };
        cold.raw_value = raw_value;
        cold.weight = total_weight - remaining_weight;
        cold.count = (data_index + 1) as u32;

        let result = raw_value.max(lower_bound).min(upper_bound);
        self.set_value(idx, row_idx, result);
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
            let idx = Self::compute_index((r, observation.x), self.len());
            let bounds = if r < observation.x {
                self.propagate_bounds_with_row::<K>(idx, r, observation.x)
            } else {
                // Diagonal (singleton) cell: no subdivisions, so no clip bounds. The
                // (NaN, NaN) sentinel makes `update_value`'s `lower >= upper` collapse
                // check false and the final `max(NaN).min(NaN)` clip a no-op.
                (f32::NAN, f32::NAN)
            };
            self.update_value::<K, _, _>(
                data_index,
                idx,
                bounds,
                r,
                observation.x,
                input,
                epsilon,
                marker_max,
            );
        }
    }

    /// Propagate bounds at `(r, s)`: reduce over the row operand `values[(r, r..s)]` (read
    /// contiguously from the row-major mirror) paired with the column operand
    /// `values[(r+1..=s, s)]` (contiguous in the column-major triangle).
    ///
    /// `idx` must be `compute_index((r, s))`. The bounds are only ever consumed by the
    /// `update_value` that immediately follows (kept in registers), so they are returned
    /// rather than stored — no per-cell bounds array is maintained.
    fn propagate_bounds_with_row<K: Kernel>(&self, idx: usize, r: usize, s: usize) -> Bounds {
        debug_assert!(r < s);
        debug_assert!(s < self.len());

        let len = s - r;
        let col_s_base = idx - r;
        debug_assert_eq!(col_s_base, s * (s + 1) / 2);
        let row_start = Self::compute_row_index((r, r), self.len());
        debug_assert!(row_start + len <= self.values_row.len());
        debug_assert!(col_s_base + s < self.values.len());
        // SAFETY: both callers guarantee `r < s < self.len()` (cell coordinates come from
        // partition indices bounded by `n_covariate`). The row operand `(r, r..s)` then ends
        // at the row-mirror index of `(r, s)`, inside row `r`'s segment, and the column
        // operand `(r+1..=s, s)` ends at the triangle index of `(s, s)`, the last entry of
        // column `s` — both in bounds, and `r < s` keeps both ranges well-formed.
        let (row, col) = unsafe {
            (
                self.values_row.get_unchecked(row_start..row_start + len),
                self.values
                    .get_unchecked(col_s_base + r + 1..col_s_base + s + 1),
            )
        };

        // Kernel-mode choice, a single comparison against the fit-level effective
        // threshold (see `Estimates::checkpoint_min_len`). Long reductions on fits with
        // real censoring run with the collapse checkpoints: one checkpoint reduce is
        // noise next to their straight-line work, and a collapsed long cell (the
        // dominant visit on weakly-dependent continuous censored fits) exits with most
        // of its work skipped. Everything else runs the straight-line kernel: short and
        // mid-length reductions cannot profit (below the first checkpoint the checks
        // cannot fire, just above it the reduces on open-bounds cells cost more than
        // the exits save — measured net regressions on weakly correlated fits, whose
        // collapsed cells concentrate at these lengths), and low-censoring fits almost
        // never collapse at all. Both dispatch inputs are deterministic and identical
        // across SIMD levels: the length is a static property of the cell, the
        // threshold a pure function of the input data.
        if len >= self.checkpoint_min_len {
            K::apply::<true>(row, col)
        } else {
            K::apply::<false>(row, col)
        }
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
