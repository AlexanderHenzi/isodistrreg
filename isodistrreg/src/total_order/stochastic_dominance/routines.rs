use crate::progress::ProgressTracker;
use crate::structures::{Decreasing, Direction, Observation};
use crate::total_order::routines;
use crate::total_order::structures::{CovariateStatistic, Partition, WeightedPartition};
use bitree::BITree;
use std::cmp::Reverse;
use std::iter::repeat_n;

pub fn accelerated_pava<R: PartialOrd + Copy, S: Copy, STOP: Fn(S) -> bool, D: Direction>(
    data_index: &mut usize,
    observations: &[Observation<usize, R, S>],
    covariate_statistics: &[CovariateStatistic],
    stop_condition: STOP,
    partitions_to_store: &mut Vec<WeightedPartition>,
    cdf: &mut Vec<f64>,
    progress: &dyn ProgressTracker,
) -> (Vec<f64>, BITree<f64>, Vec<WeightedPartition>) {
    debug_assert!(
        observations.len() - *data_index > 1,
        "Need at least two observations; one to initialize, another to start the loop",
    );

    let n_covariate = covariate_statistics.len();

    // State variables maintained as we iterate over the data
    let mut consumed_share = vec![0.0; n_covariate];
    // TODO: Keeping this data structure updated has a cost that is only worth it from n=2000-3000
    let mut consumed_weight = BITree::new_zeros(n_covariate);
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
    consumed_share: &mut [f64],
    consumed_weight: &mut BITree<f64>,
    partitions: &mut Vec<WeightedPartition>,
    observations: &[Observation<usize, R, S>],
    covariate_statistics: &[CovariateStatistic],
) {
    // Get first value for lowest y and build the antitonic regression
    let obs = &observations[*data_index];
    consumed_share[obs.x] = obs.weight / covariate_statistics[obs.x].weight;
    consumed_weight.add_at(obs.x, obs.weight);

    let total_weight = covariate_statistics.last().unwrap().cumulative_weight;
    let n_covariate = covariate_statistics.len();

    // When we initialize the partitions for increasing (resp. decreasing), the thresholds are
    // decreasing (resp. increasing) isotonic regressions. D is the ordering of the threshold.
    if !D::IS_INCREASING {
        let left_weight = covariate_statistics[obs.x].cumulative_weight;
        partitions.push(Partition {
            index: obs.x + 1,
            weight: left_weight,
            value: obs.weight / left_weight,
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
                value: obs.weight / total_weight,
            });
        } else {
            let left_weight = covariate_statistics[obs.x - 1].cumulative_weight;
            let right_weight = total_weight - left_weight;
            partitions.push(Partition {
                index: n_covariate,
                weight: right_weight,
                value: obs.weight / right_weight,
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
    obs: &Observation<usize, R, S>,
    consumed_share: &mut [f64],
    consumed_weight: &mut BITree<f64>,
    partitions: &mut Vec<WeightedPartition>,
    covariate_statistics: &[CovariateStatistic],
    partitions_to_store: &mut Vec<WeightedPartition>,
) {
    let covariate = obs.x;

    consumed_share[covariate] += obs.weight / covariate_statistics[covariate].weight;
    consumed_weight.add_at(covariate, obs.weight);

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
    routines::pool_partitions_from_right_can_reindex::<Decreasing>(partitions, !D::IS_INCREASING);

    // Accelerated extension and pooling
    match D::IS_INCREASING {
        false => add_remaining(
            covariate + 1..upper,
            partitions,
            consumed_share,
            covariate_statistics,
            true,
        ),
        true => add_remaining(
            (lower..covariate).rev(),
            partitions,
            consumed_share,
            covariate_statistics,
            false,
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

fn new_partition(
    lower: usize,
    upper_inclusive: usize,
    consumed: &BITree<f64>,
    covariate_statistics: &[CovariateStatistic],
) -> Partition<f64, f64> {
    let total_weight = if lower > 0 {
        covariate_statistics[upper_inclusive].cumulative_weight
            - covariate_statistics[lower - 1].cumulative_weight
    } else {
        covariate_statistics[upper_inclusive].cumulative_weight
    };
    let consumed_weight = range_sum_inclusive(consumed, lower, upper_inclusive);
    Partition {
        index: upper_inclusive + 1,
        weight: total_weight,
        value: consumed_weight / total_weight,
    }
}

/// Inclusive-range prefix-sum query on the BIT: sum over covariate indices `[a, b]`.
/// `ftree::BITree::prefix_sum(k, init)` returns `init + Σ_{i<k} values[i]`, i.e. strict
/// prefix, so the inclusive range `[a, b]` corresponds to `prefix_sum(b+1) − prefix_sum(a)`.
#[inline]
fn range_sum_inclusive(tree: &BITree<f64>, a: usize, b: usize) -> f64 {
    debug_assert!(a <= b);
    tree.prefix_sum(b + 1) - tree.prefix_sum(a)
}

fn add_remaining(
    range_to_add: impl Iterator<Item = usize>,
    partitions: &mut Vec<Partition<f64, f64>>,
    consumed_share: &mut [f64],
    covariate_statistics: &[CovariateStatistic],
    reindex_on_pool: bool,
) {
    // TODO: Can we track how blocks were merged to avoid having to add singletons, and instead
    //  merge larger sub blocks? (not in the worst case, but in the typical case, maybe?)
    for i in range_to_add {
        partitions.push(Partition {
            index: i + 1,
            weight: covariate_statistics[i].weight,
            value: consumed_share[i],
        });
        routines::pool_partitions_from_right_can_reindex::<Decreasing>(partitions, reindex_on_pool);
    }
}

/// Store the latest antitonic regression as represented by the partitions in the CDF.
///
/// We want to store a dense output (that we transpose later) because it is fast for predicting a
/// conditional cdf, but this is the only O(n^2) part of the algorithm so it start to dominate for
/// large inputs.
pub fn store_in_cdf<W, D: Direction>(partitions: &[Partition<W, f64>], cdf: &mut Vec<f64>) {
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
