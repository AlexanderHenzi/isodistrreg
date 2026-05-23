use crate::structures::Direction;
use crate::total_order::structures::{CovariateStatistic, Partition};
use crate::total_order::tonic_regression;
use crate::{Float, Observation};

/// Pool by weighted averaging until the partitions values follow the specified direction.
pub fn pool_partitions_from_right<D: Direction, F: Float>(parts: &mut Vec<Partition<F, F>>) {
    pool_partitions_from_right_can_reindex::<D, F>(parts, true);
}

/// Re-indexing is not necessary when the partitions are sorted in reverse order, like when
/// maintaining a partition for a decreasing (S-)IDR (an increasing set of thresholds).
pub fn pool_partitions_from_right_can_reindex<D: Direction, F: Float>(
    parts: &mut Vec<Partition<F, F>>,
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
        absorbs.weight = absorbs.weight + gets_absorbed.weight;
        absorbs.value = absorbs.value / absorbs.weight;
    }
}

/// Fast path for a single-threshold fit: a binary isotonic regression where each covariate's
/// "response" is its share of uncensored weight.
pub fn single_response<D: Direction, Y>(
    observations: Vec<Observation<usize, Y, bool, f32>>,
    covariate_statistics: Vec<CovariateStatistic>,
) -> Vec<f32> {
    let n_covariate = covariate_statistics.len();
    // observations may not be sorted by covariate
    let mut uncensored_per_covariate = vec![0.0f32; n_covariate];
    for o in observations {
        if o.observed {
            uncensored_per_covariate[o.x] += o.weight;
        }
    }
    let shares = uncensored_per_covariate
        .into_iter()
        .zip(covariate_statistics.iter())
        .map(|(uncensored, cs)| (uncensored / cs.weight, cs.weight));
    tonic_regression::algorithm_pre_sorted_deduplicated::<D::REVERSE, f32, _>(shares).collect()
}
