use crate::Observation;
use crate::structures::Direction;
use crate::total_order::structures::{CovariateStatistic, WeightedPartition};
use crate::total_order::tonic_regression::algorithm_pre_sorted_deduplicated as tonic_regression_pre_sorted;

/// Pool by weighted averaging until the partitions values follow the specified direction.
pub fn pool_partitions_from_right<D: Direction>(parts: &mut Vec<WeightedPartition>) {
    pool_partitions_from_right_can_reindex::<D>(parts, true);
}

/// Re-indexing is not necessary when the partitions are sorted in reverse order, like when
/// maintaining a partition for a decreasing (S-)IDR (an increasing set of thresholds).
pub fn pool_partitions_from_right_can_reindex<D: Direction>(
    parts: &mut Vec<WeightedPartition>,
    with_reindex: bool,
) {
    while matches!(
        parts.as_slice(),
        [.., prev, last] if prev.value.partial_cmp(&last.value).unwrap() == D::FORBIDDEN_ORDERING,
    ) {
        let gets_absorbed = parts.pop().unwrap();
        let absorbs = parts.last_mut().unwrap();

        if with_reindex {
            absorbs.index = gets_absorbed.index;
        }
        absorbs.value = absorbs.weight * absorbs.value + gets_absorbed.weight * gets_absorbed.value;
        absorbs.weight += gets_absorbed.weight;
        absorbs.value /= absorbs.weight;
    }
}

pub fn single_response<D: Direction, Y>(
    observations: Vec<Observation<usize, Y, bool>>,
    covariate_statistics: Vec<CovariateStatistic>,
) -> Vec<f64> {
    let n_covariate = covariate_statistics.len();
    // observations may not be sorted by covariate
    let mut uncensored_per_covariate = vec![0.0; n_covariate];
    for o in observations {
        if o.observed {
            uncensored_per_covariate[o.x] += o.weight;
        }
    }
    let share_uncensored_and_weight = (0..n_covariate).map(|i| {
        let total_weight = covariate_statistics[i].weight;
        (uncensored_per_covariate[i] / total_weight, total_weight)
    });
    tonic_regression_pre_sorted::<D::REVERSE>(share_uncensored_and_weight).collect()
}
