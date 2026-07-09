use crate::progress::ProgressTracker;
use crate::structures::{Decreasing, Direction, Observation};
use crate::total_order::routines;
use crate::total_order::structures::{CovariateStatistic, Partition, WeightedPartition};
use bitree::BITree;
use std::cmp::Reverse;
use std::iter::repeat_n;

/// Slack for the debug-only monotonicity asserts on emitted CDF rows and partition
/// values. Exact monotonicity is a property of the estimator in exact arithmetic, but not
/// of the f32 evaluation: batched pooling folds along shorter merge paths, and the
/// censored collapse checkpoints return midpoints of per-visit checkpoint prefixes that
/// are mutually unordered — both wobble results by rounding noise (measured elsewhere at
/// <= ~4e-7). The slack separates that noise from logic bugs, which violate monotonicity
/// by orders of magnitude more.
#[cfg(debug_assertions)]
pub(crate) const MONOTONICITY_NOISE_SLACK: f32 = 1e-5;

/// Debug-only: compare the last two stored CDF rows — the distribution function at each
/// covariate must be non-decreasing across thresholds (up to
/// [`MONOTONICITY_NOISE_SLACK`]). Downstream consumers (quantile extraction, the
/// proper-CDF gates) rely on the emitted rows forming valid CDFs.
#[cfg(debug_assertions)]
pub(crate) fn debug_assert_last_rows_monotone(cdf: &[f32], n_covariate: usize) {
    if cdf.len() < 2 * n_covariate {
        return;
    }
    let tail = &cdf[cdf.len() - 2 * n_covariate..];
    let (prev, curr) = tail.split_at(n_covariate);
    for (i, (p, q)) in prev.iter().zip(curr).enumerate() {
        assert!(
            *q >= *p - MONOTONICITY_NOISE_SLACK,
            "CDF decreased across thresholds at covariate {i}: {p} -> {q}",
        );
    }
}

/// Exact bookkeeping of the observations the classical PAVA has consumed so far.
///
/// All weight arithmetic runs in f32 (see the weight contract on `fit()`), so quotients of
/// consumed and total weight drift by ulps and cannot be trusted to land on exact 0.0/1.0.
/// The integer observation counts are exact: a covariate range whose consumed count equals
/// its total count is fully consumed, and its share is mathematically `total/total` — the
/// same combinatorial-pinning idea as the censored `CompletionIndex` and the hazard-rate
/// kernels' `remaining_observations`.
pub struct ConsumedMass {
    /// Consumed f32 weight per covariate, range-queryable.
    pub weight: BITree<f32>,
    /// Consumed observation count per covariate, range-queryable.
    count: BITree<u32>,
    /// `cumulative_count[i]` = number of observations with covariate index `<= i`.
    cumulative_count: Vec<u32>,
}

impl ConsumedMass {
    pub fn new<R, S>(observations: &[Observation<usize, R, S, f32>], n_covariate: usize) -> Self {
        let mut cumulative_count = vec![0u32; n_covariate];
        for observation in observations {
            cumulative_count[observation.x] += 1;
        }
        let mut accumulated = 0;
        for count in cumulative_count.iter_mut() {
            accumulated += *count;
            *count = accumulated;
        }
        Self {
            weight: BITree::new_zeros(n_covariate),
            count: BITree::new_zeros(n_covariate),
            cumulative_count,
        }
    }

    fn consume(&mut self, x: usize, weight: f32) {
        self.weight.add_at(x, weight);
        self.count.add_at(x, 1);
    }

    /// Has every observation with covariate index in `[lower, upper_inclusive]` been
    /// consumed? Pure integer arithmetic — exact.
    fn range_complete(&self, lower: usize, upper_inclusive: usize) -> bool {
        let total = self.cumulative_count[upper_inclusive]
            - if lower > 0 {
                self.cumulative_count[lower - 1]
            } else {
                0
            };
        self.count.prefix_sum(upper_inclusive + 1) - self.count.prefix_sum(lower) == total
    }

    /// Consumed f32 weight over covariate indices `[lower, upper_inclusive]`.
    fn range_weight(&self, lower: usize, upper_inclusive: usize) -> f32 {
        self.weight.prefix_sum(upper_inclusive + 1) - self.weight.prefix_sum(lower)
    }
}

pub fn accelerated_pava<R: PartialOrd + Copy, S: Copy, STOP: Fn(S) -> bool, D: Direction>(
    data_index: &mut usize,
    observations: &[Observation<usize, R, S, f32>],
    covariate_statistics: &[CovariateStatistic],
    stop_condition: STOP,
    partitions_to_store: &mut Vec<WeightedPartition>,
    cdf: &mut Vec<f32>,
    progress: &dyn ProgressTracker,
) -> (Vec<f32>, ConsumedMass, Vec<WeightedPartition>) {
    debug_assert!(
        observations.len() - *data_index > 1,
        "Need at least two observations; one to initialize, another to start the loop",
    );
    // The sort contract the run gathering below silently depends on: observations sorted
    // by response, and within one response the non-stop observations (events) form one
    // contiguous prefix with strictly ascending covariates (duplicates aggregated by
    // preprocessing). Without contiguity a run would end early and the batched update
    // would evaluate a different estimator than the per-arrival definition.
    #[cfg(debug_assertions)]
    for w in observations[*data_index..].windows(2) {
        let (a, b) = (&w[0], &w[1]);
        debug_assert!(a.y <= b.y, "observations must be sorted by response");
        if a.y == b.y {
            debug_assert!(
                !stop_condition(a.observed) || stop_condition(b.observed),
                "a threshold's non-stop observations must precede its stop observations",
            );
            if !stop_condition(a.observed) && !stop_condition(b.observed) {
                debug_assert!(
                    a.x < b.x,
                    "a threshold's events must be strictly ascending in covariate",
                );
            }
        }
    }

    let n_covariate = covariate_statistics.len();

    // State variables maintained as we iterate over the data — all in f32, matching the
    // post-preprocessing convention used by the total-order context structs, plus the exact
    // consumed-observation counts that pin fully-consumed shares to exactly 1.0.
    let mut consumed_share = vec![0.0; n_covariate];
    // TODO: Keeping this data structure updated has a cost that is only worth it from n=2000-3000
    let mut consumed_weight = ConsumedMass::new(observations, n_covariate);
    // TODO: Should this be a linked list or other data structure that makes local modifications
    //  cheap? It should be very fast to scan, though...
    let mut partitions = Vec::with_capacity(n_covariate);

    initialize_partitions::<_, _, D::REVERSE>(
        data_index,
        &mut consumed_share,
        &mut consumed_weight,
        &mut partitions,
        observations,
        covariate_statistics,
    );
    let mut active_threshold = observations[*data_index].y;
    *data_index += 1;

    let last_threshold = observations.last().unwrap().y;
    // Each iteration completes one threshold: its events form one contiguous
    // ascending-covariate run (preprocessing sorts each threshold by covariate,
    // aggregates duplicates, and puts a threshold's events before its stop
    // observations), applied as a single batched update, followed by the threshold's
    // CDF row. Only the very first iteration can instead begin mid-threshold:
    // `initialize_partitions` consumed the first observation.
    loop {
        if observations[*data_index].y == active_threshold
            && !stop_condition(observations[*data_index].observed)
        {
            let mut run_end = *data_index + 1;
            while run_end < observations.len() {
                let next = &observations[run_end];
                if next.y != active_threshold || stop_condition(next.observed) {
                    break;
                }
                run_end += 1;
            }
            update_threshold::<_, _, D::REVERSE>(
                &observations[*data_index..run_end],
                &mut consumed_share,
                &mut consumed_weight,
                &mut partitions,
                covariate_statistics,
                partitions_to_store,
            );
            *data_index = run_end;
        }

        // A stop observation inside the current threshold ends the prefix before the
        // threshold's row is stored.
        if observations[*data_index].y == active_threshold {
            debug_assert!(stop_condition(observations[*data_index].observed));
            return (consumed_share, consumed_weight, partitions);
        }

        #[cfg(debug_assertions)]
        let cdf_len_before = cdf.len();
        store_in_cdf::<_, D>(&partitions, cdf);
        #[cfg(debug_assertions)]
        {
            // The partition stack must tile the covariate grid exactly — one full row
            // per threshold — and successive rows must form valid per-covariate CDFs.
            debug_assert_eq!(
                cdf.len() - cdf_len_before,
                n_covariate,
                "a threshold's partitions did not tile the covariate grid",
            );
            debug_assert_last_rows_monotone(cdf, n_covariate);
        }
        progress.increment();
        active_threshold = observations[*data_index].y;

        if active_threshold == last_threshold {
            break;
        }
    }

    (consumed_share, consumed_weight, partitions)
}

/// Set up the first partitioning.
///
/// Based on just a single observation (lowest threshold, one of the (potentially multiple
/// non-unique response values)) we build the simple partitioning of just one or two partitions.
/// This special case initializes the partitions and speeds up the first iteration.
///
/// If the IDR needs to be stochastically increasing (nondecreasing) w.r.t. the covariate, each
/// threshold is a decreasing (nonincreasing) regression. If the (S-)IDR is increasing, this method
/// should be called with D = Decreasing.
fn initialize_partitions<R, S, D: Direction>(
    data_index: &mut usize,
    consumed_share: &mut [f32],
    consumed_weight: &mut ConsumedMass,
    partitions: &mut Vec<WeightedPartition>,
    observations: &[Observation<usize, R, S, f32>],
    covariate_statistics: &[CovariateStatistic],
) {
    // Get first value for lowest y and build the antitonic regression
    let obs = &observations[*data_index];
    consumed_weight.consume(obs.x, obs.weight);
    // The share is exactly 1.0 when this is the covariate's only observation (the f32
    // quotient of separately accumulated sums can sit ulps off).
    consumed_share[obs.x] = if consumed_weight.range_complete(obs.x, obs.x) {
        1.0
    } else {
        obs.weight / covariate_statistics[obs.x].weight
    };

    let total_weight = covariate_statistics.last().unwrap().cumulative_weight;
    let n_covariate = covariate_statistics.len();

    // The first block's value is its consumed share; pinned to exactly 1.0 when the block
    // holds no other observation, like `new_partition`.
    let block_share = |lower: usize, upper_inclusive: usize, block_weight: f32| {
        if consumed_weight.range_complete(lower, upper_inclusive) {
            1.0
        } else {
            obs.weight / block_weight
        }
    };

    // When we initialize the partitions for increasing (resp. decreasing), the thresholds are
    // decreasing (resp. increasing) isotonic regressions. D is the ordering of the threshold.
    if !D::IS_INCREASING {
        let left_weight = covariate_statistics[obs.x].cumulative_weight;
        partitions.push(Partition {
            index: obs.x + 1,
            weight: left_weight,
            value: block_share(0, obs.x, left_weight),
        });
        // The first antitonic regression may have a second partition, if the value isn't the last
        if obs.x < n_covariate - 1 {
            partitions.push(Partition {
                index: n_covariate,
                weight: total_weight - left_weight,
                value: 0.0,
            });
        }
    } else {
        if obs.x == 0 {
            partitions.push(Partition {
                index: n_covariate,
                weight: total_weight,
                value: block_share(0, n_covariate - 1, total_weight),
            });
        } else {
            let left_weight = covariate_statistics[obs.x - 1].cumulative_weight;
            let right_weight = total_weight - left_weight;
            partitions.push(Partition {
                index: n_covariate,
                weight: right_weight,
                value: block_share(obs.x, n_covariate - 1, right_weight),
            });
            partitions.push(Partition {
                index: obs.x,
                weight: left_weight,
                value: 0.0,
            });
        }
    }
    debug_assert!(
        partitions.iter().is_sorted_by_key(|p| Reverse(p.value)),
        "For increasing regressions, we store the partitions in reverse",
    );
}

/// Apply one threshold's event run — `run`, tied on the response, covariates strictly
/// ascending by the preprocessing sort and aggregation — as one batched PAVA update. A
/// run of one is exactly the classical per-arrival update step — split-inheritance,
/// pooling, and the trailing flush as the accelerated extension — so there is no
/// separate single-arrival path.
///
/// Per-arrival processing re-adds a split block's deferred remainder once per arrival,
/// and the next tied arrival splits the freshly rebuilt structure again — a run of k
/// tied observations re-adds overlapping remainders up to k times (measured 10-33x the
/// covariate-axis width per threshold on weakly-correlated tied-response inputs, where
/// the update loop is >= 92% of the fit). Only the state after the whole run is ever
/// read (`store_in_cdf` runs once per threshold), and the classical PAVA solution is
/// insertion-order independent, so the run is applied in a single pass that re-adds
/// every covariate of the affected range at most once:
///
/// - The consumed bookkeeping for the whole run happens first, so every partition
///   value computed below reads threshold-end state. This is equivalent to the
///   per-arrival flow: each value it computes spans covariates on the processed side
///   of its arrival, and the later arrivals of the run lie strictly beyond them.
/// - Arrivals are processed in storage order (see below), so a block built during the
///   run spans covariates on the processed side of the arrivals and no later arrival
///   of the run splits it again.
/// - A split block's remainder is re-added only when the next arrival must pass it
///   (or by the trailing flush). An arrival inside such a remainder enters through the
///   same `add_remaining` stack pass as its remainder neighbors — its consumed share
///   is already final, so this is the canonical classical-PAVA insertion.
/// - An arrival inside an untouched original block splits it by inheritance, the
///   per-arrival move. Untouched blocks between arrivals are restored verbatim without
///   pooling, like the per-arrival flow's blind restore: their values and their
///   boundary's validity are unchanged.
///
/// The two update directions are mirror images of each other, and the reversed
/// partition storage of the `Increasing` direction absorbs the mirroring: in both
/// directions the arrivals are processed in STORAGE order (ascending covariate for
/// `Decreasing`, descending for `Increasing` — always from the stack's live end into
/// the drained store), a split block keeps its storage-near part, and the pending
/// region lies between the stack top and the next stored block. The body below is
/// written in those storage terms, tracking the pending region by its `(near, far)`
/// covariate edges; the only direction-dependent lines are the coordinate
/// conversions, each folded at compile time.
///
/// Pooled block means accumulate along a different (shorter) merge path than the
/// per-arrival flow, so outputs can differ within rounding noise — the same contract
/// as `add_remaining`'s equal-value coalescing, gated by the f64-reference
/// differential tests.
pub fn update_threshold<R, S, D: Direction>(
    run: &[Observation<usize, R, S, f32>],
    consumed_share: &mut [f32],
    consumed_weight: &mut ConsumedMass,
    partitions: &mut Vec<WeightedPartition>,
    covariate_statistics: &[CovariateStatistic],
    partitions_to_store: &mut Vec<WeightedPartition>,
) {
    debug_assert!(!run.is_empty());
    debug_assert!(run.windows(2).all(|w| w[0].x < w[1].x));

    // Threshold-end bookkeeping for the whole run, with the exact-count pinning that
    // keeps mathematically flat CDF rows exactly at 1.0.
    for observation in run {
        let covariate = observation.x;
        consumed_weight.consume(covariate, observation.weight);
        consumed_share[covariate] = if consumed_weight.range_complete(covariate, covariate) {
            1.0
        } else {
            consumed_share[covariate] + observation.weight / covariate_statistics[covariate].weight
        };
    }

    // Re-add covariates `[lo, hi)` as singleton blocks, in storage order.
    fn re_add<D: Direction>(
        lo: usize,
        hi: usize,
        partitions: &mut Vec<WeightedPartition>,
        consumed_share: &[f32],
        covariate_statistics: &[CovariateStatistic],
    ) {
        if !D::IS_INCREASING {
            add_remaining::<true, _>(lo..hi, partitions, consumed_share, covariate_statistics);
        } else {
            add_remaining::<false, _>(
                (lo..hi).rev(),
                partitions,
                consumed_share,
                covariate_statistics,
            );
        }
    }

    // First arrival in storage order: the per-arrival split-inheritance move on its
    // containing block, keeping the block's storage-near part.
    let first_x = if !D::IS_INCREASING {
        run[0].x
    } else {
        run[run.len() - 1].x
    };
    let (partition_index, (lower, upper)) = find_partition_bounds::<_, _, D>(first_x, partitions);
    debug_assert!(lower <= first_x && first_x < upper);
    partitions_to_store.extend(partitions.drain(partition_index + 1..));
    partitions[partition_index] = if !D::IS_INCREASING {
        new_partition(lower, first_x, consumed_weight, covariate_statistics)
    } else {
        new_partition(first_x, upper - 1, consumed_weight, covariate_statistics)
    };
    routines::pool_partitions_from_right_can_reindex::<Decreasing, _>(
        partitions,
        !D::IS_INCREASING,
    );

    // Covariates strictly between `near` and `far` (in storage order: `near` abuts the
    // stack top, `far` the next stored block) are not on the stack yet: the deferred
    // remainder of the last split block when `pending_is_remainder`, otherwise an
    // untouched original block. `partitions_to_store[store_cursor..]` holds the blocks
    // not yet reached, in storage order.
    let (mut near, mut far) = if !D::IS_INCREASING {
        (first_x + 1, upper)
    } else {
        (first_x, lower)
    };
    let mut pending_is_remainder = true;
    let mut store_cursor = 0usize;

    for i in 1..run.len() {
        let observation = if !D::IS_INCREASING {
            &run[i]
        } else {
            &run[run.len() - 1 - i]
        };
        let x = observation.x;

        // Advance the stack up to the pending region containing x.
        while if !D::IS_INCREASING { far <= x } else { x < far } {
            if pending_is_remainder {
                let (lo, hi) = if !D::IS_INCREASING {
                    (near, far)
                } else {
                    (far, near)
                };
                re_add::<D>(lo, hi, partitions, consumed_share, covariate_statistics);
            } else {
                // The pending region is the untouched block just pulled from the store.
                partitions.push(partitions_to_store[store_cursor - 1].clone());
            }
            near = far;
            far = if !D::IS_INCREASING {
                partitions_to_store[store_cursor].index
            } else {
                debug_assert_eq!(partitions_to_store[store_cursor].index, near);
                partitions_to_store
                    .get(store_cursor + 1)
                    .map_or(0, |block| block.index)
            };
            pending_is_remainder = false;
            store_cursor += 1;
        }

        // The storage-order position just past the arrival; also the pending region's
        // near edge once the arrival is on the stack.
        let past_x = if !D::IS_INCREASING { x + 1 } else { x };

        if pending_is_remainder {
            // x lies in a deferred remainder: one storage-order stack pass re-adds the
            // covariates between the stack top and x, and inserts x itself.
            let (lo, hi) = if !D::IS_INCREASING {
                (near, past_x)
            } else {
                (past_x, near)
            };
            re_add::<D>(lo, hi, partitions, consumed_share, covariate_statistics);
        } else {
            // x lies in an untouched original block: split-inheritance, exactly like
            // the first arrival.
            partitions.push(if !D::IS_INCREASING {
                new_partition(near, x, consumed_weight, covariate_statistics)
            } else {
                new_partition(x, near - 1, consumed_weight, covariate_statistics)
            });
            routines::pool_partitions_from_right_can_reindex::<Decreasing, _>(
                partitions,
                !D::IS_INCREASING,
            );
        }
        near = past_x;
        pending_is_remainder = true;
    }

    // Trailing remainder of the last split block, then the blind restore of the
    // untouched far side — the per-arrival flow's closing moves.
    let (lo, hi) = if !D::IS_INCREASING {
        (near, far)
    } else {
        (far, near)
    };
    re_add::<D>(lo, hi, partitions, consumed_share, covariate_statistics);
    partitions.extend(partitions_to_store.drain(store_cursor..));
    partitions_to_store.clear();

    debug_assert!(partitions.is_sorted_by_key(|p| Reverse(p.value)));
    #[cfg(debug_assertions)]
    {
        // The restored stack must tile the covariate grid in storage order (ascending
        // block boundaries for the decreasing direction, descending for the increasing
        // direction's reversed storage) — a gap or overlap here silently corrupts every
        // later threshold. Pooled means divide by the summed block weight and the pool
        // trichotomy reads finite shares, so both must stay valid through the batched
        // splits/re-adds/restores.
        if !D::IS_INCREASING {
            debug_assert!(partitions.windows(2).all(|w| w[0].index < w[1].index));
            debug_assert_eq!(partitions.last().unwrap().index, covariate_statistics.len(),);
        } else {
            debug_assert!(partitions.windows(2).all(|w| w[0].index > w[1].index));
            debug_assert_eq!(
                partitions.first().unwrap().index,
                covariate_statistics.len(),
            );
        }
        debug_assert!(
            partitions
                .iter()
                .all(|p| p.weight > 0.0 && p.value.is_finite() && p.value >= 0.0),
            "partition weights must stay positive and values finite non-negative shares",
        );
    }
}

pub fn find_partition_bounds<W, V, D: Direction>(
    covariate_index: usize,
    partitions: &[Partition<W, V>],
) -> (usize, (usize, usize)) {
    let (partition_index, lower) = if !D::IS_INCREASING {
        debug_assert!(partitions.is_sorted_by_key(|p| p.index));

        // Find partition to which covariate value belongs
        let found = partitions.binary_search_by_key(&covariate_index, |p| p.index);
        let partition_index = match found {
            Ok(index) => index + 1,
            Err(index) => index,
        };
        // Inclusive
        let lower = if partition_index > 0 {
            partitions[partition_index - 1].index
        } else {
            debug_assert_eq!(partition_index, 0);
            0
        };

        (partition_index, lower)
    } else {
        debug_assert!(partitions.is_sorted_by_key(|p| Reverse(p.index)));

        // Find partition to which covariate value belongs
        let found =
            partitions.binary_search_by_key(&Reverse(covariate_index), |p| Reverse(p.index));
        let partition_index = found.unwrap_or_else(|i| i) - 1;
        // Inclusive
        let lower = if partition_index < partitions.len() - 1 {
            partitions[partition_index + 1].index
        } else {
            debug_assert_eq!(partition_index, partitions.len() - 1);
            0
        };

        (partition_index, lower)
    };

    // Exclusive
    let upper = partitions[partition_index].index;

    (partition_index, (lower, upper))
}

/// Build the partition covering covariates `[lower, upper_inclusive]`, deriving its total
/// weight from the cumulative weights and its value from the consumed-mass BIT. A fully
/// consumed block (exact integer-count check) is pinned to exactly 1.0: the value is
/// mathematically `total/total`, while the two independently accumulated f32 sums can each
/// drift by ulps.
fn new_partition(
    lower: usize,
    upper_inclusive: usize,
    consumed: &ConsumedMass,
    covariate_statistics: &[CovariateStatistic],
) -> Partition<f32, f32> {
    let total_weight = if lower > 0 {
        covariate_statistics[upper_inclusive].cumulative_weight
            - covariate_statistics[lower - 1].cumulative_weight
    } else {
        covariate_statistics[upper_inclusive].cumulative_weight
    };
    let value = if consumed.range_complete(lower, upper_inclusive) {
        1.0
    } else {
        consumed.range_weight(lower, upper_inclusive) / total_weight
    };
    Partition {
        index: upper_inclusive + 1,
        weight: total_weight,
        value,
    }
}

/// Re-add the covariates in `range_to_add` as singleton blocks, pooling each into the
/// running partition stack so it stays a `Decreasing`-by-value antitonic regression.
///
/// This is the accelerated PAVA's hottest loop on weakly-correlated data (a split block
/// re-adds its whole interior). Rather than `push` a singleton and immediately call the
/// generic pool — which pops it straight back out whenever it merges — the merge is
/// accumulated in registers and the block is pushed once. `REINDEX` (compile-time,
/// mirroring the generic pool's runtime flag) selects which merged endpoint's `index`
/// (block boundary) survives.
///
/// Besides the strict `Decreasing` violation (pooled by weighted mean, operand order
/// matching the generic pool), an *equal*-valued neighbour is coalesced: it carries the
/// same value, so the merged block keeps it bit-for-bit and only sums the weight. This is
/// the important part — the regression is frequently a step function whose runs a strict
/// pool leaves as thousands of singletons (a long consumed/unconsumed plateau); folding
/// them keeps the list near its true block count, shrinking every later suffix copy, store
/// and scan. The emitted CDF is the same step function; only later pooled means round
/// differently (a coalesced weight vs. a chain of pair merges), all within the f32 noise
/// floor — gated by the f64-reference differential tests. Values are finite shares in
/// `[0, 1]`, so `<`/`==` match the `partial_cmp` trichotomy the generic pool used.
fn add_remaining<const REINDEX: bool, I: Iterator<Item = usize>>(
    range_to_add: I,
    partitions: &mut Vec<WeightedPartition>,
    consumed_share: &[f32],
    covariate_statistics: &[CovariateStatistic],
) {
    for i in range_to_add {
        let mut acc = Partition {
            index: i + 1,
            weight: covariate_statistics[i].weight,
            value: consumed_share[i],
        };
        while let Some(prev) = partitions.last() {
            if prev.value < acc.value {
                let prev = partitions.pop().unwrap();
                let pooled_weight = prev.weight + acc.weight;
                acc.value = (prev.weight * prev.value + acc.weight * acc.value) / pooled_weight;
                acc.weight = pooled_weight;
                if !REINDEX {
                    acc.index = prev.index;
                }
            } else if prev.value == acc.value {
                let prev = partitions.pop().unwrap();
                acc.weight += prev.weight;
                if !REINDEX {
                    acc.index = prev.index;
                }
            } else {
                break;
            }
        }
        partitions.push(acc);
    }
}

/// Store the latest antitonic regression as represented by the partitions in the CDF.
///
/// All quantities are f32 throughout; we still clamp before saving because we sometimes shave
/// off a tiny epsilon above 1.0.
pub fn store_in_cdf<W, D: Direction>(partitions: &[Partition<W, f32>], cdf: &mut Vec<f32>) {
    let stored_in_reverse = !D::IS_INCREASING;
    if !stored_in_reverse {
        let first = partitions.first().unwrap();
        cdf.extend(repeat_n(first.value.clamp(0.0, 1.0), first.index));
        for l in 1..partitions.len() {
            let partition_len = partitions[l].index - partitions[l - 1].index;
            cdf.extend(repeat_n(partitions[l].value.clamp(0.0, 1.0), partition_len));
        }
    } else {
        let first = partitions.last().unwrap();
        cdf.extend(repeat_n(first.value.clamp(0.0, 1.0), first.index));
        for l in (0..partitions.len() - 1).rev() {
            let partition_len = partitions[l].index - partitions[l + 1].index;
            cdf.extend(repeat_n(partitions[l].value.clamp(0.0, 1.0), partition_len));
        }
    }
}
