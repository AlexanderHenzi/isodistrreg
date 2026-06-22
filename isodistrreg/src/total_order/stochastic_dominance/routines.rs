use crate::progress::ProgressTracker;
use crate::structures::{Decreasing, Direction, Observation};
use crate::total_order::routines;
use crate::total_order::structures::{CovariateStatistic, Partition, WeightedPartition};
use bitree::BITree;
use std::cmp::Reverse;
use std::iter::repeat_n;

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
    loop {
        while observations[*data_index].y == active_threshold {
            let observation = &observations[*data_index];

            // TODO: We could do a single, compacted pava update step in which we mark all
            //  partition elements in which a new value lands as "dirty" and don't end up splitting
            //  some partition elements twice.

            if stop_condition(observation.observed) {
                return (consumed_share, consumed_weight, partitions);
            }

            classical_pava_update_step::<_, _, D::REVERSE>(
                observation,
                &mut consumed_share,
                &mut consumed_weight,
                &mut partitions,
                covariate_statistics,
                partitions_to_store,
            );

            *data_index += 1;
        }

        store_in_cdf::<_, D>(&partitions, cdf);
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

/// If the IDR is increasing, should be called with D = Decreasing.
pub fn classical_pava_update_step<R, S, D: Direction>(
    obs: &Observation<usize, R, S, f32>,
    consumed_share: &mut [f32],
    consumed_weight: &mut ConsumedMass,
    partitions: &mut Vec<WeightedPartition>,
    covariate_statistics: &[CovariateStatistic],
    partitions_to_store: &mut Vec<WeightedPartition>,
) {
    let covariate = obs.x;

    consumed_weight.consume(covariate, obs.weight);
    // Once every observation of the covariate is consumed, the share is mathematically
    // `total/total` — exactly 1.0; the f32 running sum can sit ulps off, and rows that are
    // mathematically flat at 1.0 must not wobble (`Fit::assert_consistent` requires every
    // per-covariate CDF row to be sorted along the thresholds).
    consumed_share[covariate] = if consumed_weight.range_complete(covariate, covariate) {
        1.0
    } else {
        consumed_share[covariate] + obs.weight / covariate_statistics[covariate].weight
    };

    // In which partition does the new observation fall?
    let (partition_index, (lower, upper)) = find_partition_bounds::<_, _, D>(covariate, partitions);
    debug_assert!(lower <= covariate && covariate < upper);

    // Store right part of partitions
    // TODO: Can we avoid this copy? `partitions` needs to be very fast to scan for cdf writing...
    partitions_to_store.extend(partitions.drain(partition_index + 1..));

    // Overwrite partition in which it falls with a new partition that includes the latest
    // observation
    let (new_lower, new_upper_inclusive) = match D::IS_INCREASING {
        false => (lower, covariate),
        true => (covariate, upper - 1),
    };
    partitions[partition_index] = new_partition(
        new_lower,
        new_upper_inclusive,
        consumed_weight,
        covariate_statistics,
    );

    // Pooling toward the left (of sorted or reverse sorted partitions), merging high values
    routines::pool_partitions_from_right_can_reindex::<Decreasing, _>(
        partitions,
        !D::IS_INCREASING,
    );

    // Accelerated extension and pooling
    match D::IS_INCREASING {
        false => add_remaining::<true, _>(
            covariate + 1..upper,
            partitions,
            consumed_share,
            covariate_statistics,
        ),
        true => add_remaining::<false, _>(
            (lower..covariate).rev(),
            partitions,
            consumed_share,
            covariate_statistics,
        ),
    }

    // Restore right-most partitions (low values)
    partitions.append(partitions_to_store);

    debug_assert!(partitions.is_sorted_by_key(|p| Reverse(p.value)));
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
