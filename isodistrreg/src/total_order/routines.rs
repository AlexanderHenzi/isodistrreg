use crate::structures::Direction;
use crate::total_order::structures::{CovariateStatistic, Partition};
use crate::total_order::tonic_regression;
use crate::{Float, Observation};

/// Pool by weighted averaging until the partitions values follow the specified direction.
///
/// Every merge must involve positive total weight (debug-asserted). This holds for the
/// stochastic-dominance and tonic-regression kernels: zero-weight blocks there either are
/// never pushed (`tonic_regression`) or always carry the value 0.0 (the SD PAVA guards),
/// so two of them compare `Equal` and the strict pooling condition never merges them.
/// Kernels whose block weights can reach exactly zero with differing values (the
/// hazard-rate at-risk masses) must use
/// [`pool_partitions_from_right_zero_weight_blocks`] instead.
pub fn pool_partitions_from_right<D: Direction, F: Float>(parts: &mut Vec<Partition<F, F>>) {
    pool_impl::<D, F, false>(parts, true);
}

/// Re-indexing is not necessary when the partitions are sorted in reverse order, like when
/// maintaining a partition for a decreasing (S-)IDR (an increasing set of thresholds).
///
/// Same positive-pooled-weight requirement as [`pool_partitions_from_right`].
pub fn pool_partitions_from_right_can_reindex<D: Direction, F: Float>(
    parts: &mut Vec<Partition<F, F>>,
    with_reindex: bool,
) {
    pool_impl::<D, F, false>(parts, with_reindex);
}

/// Variant of [`pool_partitions_from_right`] for kernels whose block weights can be
/// exactly zero while their values differ — the hazard-rate kernels, where the weight is
/// an at-risk mass that a censored observation can consume entirely (and, uncensored,
/// `survival² · w` can underflow). Pools a zero-total-weight pair to the midpoint of the
/// two values instead of computing 0/0 = NaN.
pub fn pool_partitions_from_right_zero_weight_blocks<D: Direction, F: Float>(
    parts: &mut Vec<Partition<F, F>>,
) {
    pool_impl::<D, F, true>(parts, true);
}

/// `ZERO_WEIGHT_BLOCKS` is a const so the fast variants compile to the plain weighted
/// average with no extra comparison in the merge loop.
#[inline]
fn pool_impl<D: Direction, F: Float, const ZERO_WEIGHT_BLOCKS: bool>(
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
        let pooled_weight = absorbs.weight + gets_absorbed.weight;
        debug_assert!(
            ZERO_WEIGHT_BLOCKS || pooled_weight > F::zero(),
            "pooled a zero-total-weight pair in a kernel that guarantees positive block weights",
        );
        absorbs.value = if !ZERO_WEIGHT_BLOCKS || pooled_weight > F::zero() {
            (absorbs.weight * absorbs.value + gets_absorbed.weight * gets_absorbed.value)
                / pooled_weight
        } else {
            // Two zero-weight blocks (e.g. covariate groups whose remaining mass was
            // fully censored away): the least-squares criterion is indifferent, any
            // value between the two satisfies the ordering — take the midpoint
            // deterministically instead of computing 0/0 = NaN.
            (absorbs.value + gets_absorbed.value) / (F::one() + F::one())
        };
        absorbs.weight = pooled_weight;
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
    let mut any_censored = vec![false; n_covariate];
    for o in observations {
        debug_assert!(o.weight > 0.0, "preprocessing drops zero-weight observations");
        if o.observed {
            uncensored_per_covariate[o.x] += o.weight;
        } else {
            any_censored[o.x] = true;
        }
    }
    let shares = uncensored_per_covariate
        .into_iter()
        .zip(any_censored)
        .zip(covariate_statistics.iter())
        .map(|((uncensored, censored), cs)| {
            let share = if !censored && uncensored > 0.0 {
                // All of this covariate's mass is observed: the share is exactly
                // total/total = 1. Summing numerator and denominator separately can
                // leave the ratio a few ulps off, which the exact proper-CDF gate in
                // `prediction::mean` must not see.
                1.0
            } else {
                uncensored / cs.weight
            };
            (share, cs.weight)
        });
    tonic_regression::algorithm_pre_sorted_deduplicated::<D::REVERSE, f32, _>(shares).collect()
}
