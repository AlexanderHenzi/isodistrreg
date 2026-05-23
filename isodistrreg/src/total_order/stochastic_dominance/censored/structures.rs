use crate::structures::Observation;
use crate::total_order;
use crate::total_order::structures::CovariateStatistic;
use bitree::BITree;

pub type Partition = total_order::structures::Partition<(), ()>;

impl Partition {
    pub fn new(index: usize) -> Self {
        Self {
            index,
            weight: (),
            value: (),
        }
    }
}

/// Preprocessed input for the censored stochastic-dominance fast algorithm.
///
/// Holds the deduplicated problem with f32 weights — preprocessing sums incoming weights
/// and downcasts at its output boundary so the entire post-preprocessing algorithm body runs
/// in f32. Covariate values stay in their input precision `X` and threshold keys in `Y`
/// (they're only used for sort/compare).
///
/// `observations[i].y` is an index into `thresholds`, not a raw response value.
#[derive(Clone, Debug, PartialEq)]
pub struct CensoredSdContext<X, Y> {
    /// Holds `n` items sorted by response from low to high, uncensored < censored (and covariate
    /// for stable tests).
    pub observations: Vec<Observation<usize, usize, bool, f32>>,
    /// Holds information about each covariate
    pub covariate_statistics: Vec<CovariateStatistic>,
    /// Unique covariate values, sorted increasing.
    pub unique_covariates: Vec<X>,
    /// Only thresholds that have at least one uncensored observation, sorted increasing.
    pub thresholds: Vec<Y>,
}

impl<X, Y> CensoredSdContext<X, Y> {
    pub fn n(&self) -> usize {
        self.observations.len()
    }
    pub fn n_covariate(&self) -> usize {
        self.covariate_statistics.len()
    }
    pub fn n_threshold(&self) -> usize {
        self.thresholds.len()
    }
}

/// The fields of a survival/Kaplan-Meier computation that the hot inner loop of `propagate_bounds`
/// does not touch. Stored in a parallel vector to `Estimates::values` so that the hot loop iterates
/// over a tightly packed `[f32]` instead of fetching a 6-field struct per access.
#[derive(Clone, Debug)]
pub struct SurvivalComputationCold {
    /// The Kaplan-Meier estimator (based on `count` samples). Updated incrementally inside
    /// `update_value`; clipped to the bounds and the result stored as the corresponding entry in
    /// `Estimates::values`.
    pub raw_value: f32,
    /// Sum of weight of included observations.
    pub weight: f32,
    /// Number of samples included in the computation of the raw value.
    pub count: usize,
    /// Lower bound based on all subdivisions of the interval - NAN for singletons.
    pub lower_bound: f32,
    /// Upper bound based on all subdivisions of the interval - NAN for singletons.
    pub upper_bound: f32,
}

/// Stores partial survival / Kaplan-Meier estimator computations for all sub intervals.
///
/// Hot/cold split: the clipped `value` of every (r, s) lives in `values` so that the hottest
/// reader (`Estimates::propagate_bounds`) walks a contiguous `[f32]` slice; the rest of the
/// per-interval state lives in `cold` at the same index.
#[derive(Debug)]
pub struct Estimates {
    len: usize,
    /// Clipped K-M value of each (r, s). Same triangular layout as `cold`. Indexed via
    /// `compute_index((r, s), len)` so `(r, i)` for `i` increasing is sequential within column `i`,
    /// and `(i+1, s)` for `i` increasing is sequential within column `s`.
    pub values: Vec<f32>,
    /// Cold half of each (r, s). Indexed identically to `values`.
    pub cold: Vec<SurvivalComputationCold>,
}

impl Estimates {
    pub fn new(n: usize, start_count: usize) -> Self {
        let total = n * (n + 1) / 2;
        let values = vec![1.0; total];
        let mut cold = vec![
            SurvivalComputationCold {
                raw_value: 1.0,
                weight: 0.0,
                count: start_count,
                lower_bound: 1.0,
                upper_bound: 1.0,
            };
            total
        ];
        // Diagonal entries don't have meaningful bounds (no subdivisions); mark with NAN so the
        // bounds-equality short-circuit in `update_value` would only fire on actually-collapsed
        // intervals.
        for i in 0..n {
            let idx = Self::compute_index((i, i), n);
            cold[idx].lower_bound = f32::NAN;
            cold[idx].upper_bound = f32::NAN;
        }
        Self {
            len: n,
            values,
            cold,
        }
    }
    /// Index into the array - `r` is inclusive, `s` is inclusive
    #[inline(always)]
    pub fn compute_index((r, s): (usize, usize), len: usize) -> usize {
        debug_assert!(r <= s);
        debug_assert!(s < len);

        s * (s + 1) / 2 + r
    }
    pub fn len(&self) -> usize {
        self.len
    }

    /// Read the clipped value at (r, s).
    #[inline(always)]
    pub fn value(&self, r: usize, s: usize) -> f32 {
        self.values[Self::compute_index((r, s), self.len)]
    }
    /// Mutable access to the cold half at (r, s).
    #[inline(always)]
    pub fn cold_mut(&mut self, r: usize, s: usize) -> &mut SurvivalComputationCold {
        &mut self.cold[Self::compute_index((r, s), self.len)]
    }
    /// Mutable access to both halves at (r, s).
    #[inline(always)]
    pub fn entry_mut(&mut self, r: usize, s: usize) -> (&mut f32, &mut SurvivalComputationCold) {
        let idx = Self::compute_index((r, s), self.len);
        (&mut self.values[idx], &mut self.cold[idx])
    }
}

impl Estimates {
    /// Bridge from the uncensored-prefix `accelerated_pava` into the generalized-PAVA hot loop.
    /// Both phases run in f32; this just unpacks the BITree's prefix sums into a plain Vec.
    pub fn from_partial_uncensored_solution(
        consumed_weight: BITree<f32>,
        covariate_statistics: &[CovariateStatistic],
        start_count: usize,
    ) -> Self {
        let len = consumed_weight.len();
        debug_assert_eq!(covariate_statistics.len(), len);

        let mut estimates = Estimates::new(len, start_count);

        let consumed_weight_plain: Vec<f32> = Vec::from(consumed_weight);

        let mut index = 0;
        for s in 0..len {
            let weight_consumed = consumed_weight_plain[s];

            // Filling in the rest of the triangle
            let (prev_values, curr_values) = estimates.values[index - s..].split_at_mut(s);
            let (prev_cold, curr_cold) = estimates.cold[index - s..].split_at_mut(s);
            for r in 0..s {
                let prev_value = prev_values[r];
                let prev_cold = &prev_cold[r];
                let curr_value = &mut curr_values[r];
                let curr_cold = &mut curr_cold[r];

                curr_cold.weight = prev_cold.weight + weight_consumed;
                let previous_full_weight = if r > 0 {
                    covariate_statistics[s - 1].cumulative_weight
                        - covariate_statistics[r - 1].cumulative_weight
                } else {
                    covariate_statistics[s - 1].cumulative_weight
                };
                let current_full_weight = covariate_statistics[s].weight;
                let share = weight_consumed / current_full_weight;
                let total_full_weight = previous_full_weight + current_full_weight;
                if total_full_weight > 0.0 {
                    curr_cold.raw_value = (previous_full_weight * prev_value
                        + current_full_weight * (1.0 - share))
                        / total_full_weight;
                } else {
                    // This covariate corresponds only has censored observations that were
                    // discarded, because they occur below any other uncensored observation, we
                    // simply replicate the value on the left instead of computing as an average
                    curr_cold.raw_value = prev_value;
                }
                debug_assert!(curr_cold.raw_value.is_finite());
                *curr_value = curr_cold.raw_value;
            }
            index += s;

            // Diagonal item last
            let share = weight_consumed / covariate_statistics[s].weight;
            let diag_cold = &mut estimates.cold[index];
            diag_cold.raw_value = 1.0 - share;
            diag_cold.weight = weight_consumed;
            estimates.values[index] = diag_cold.raw_value;
            index += 1;
        }

        estimates
    }
}
