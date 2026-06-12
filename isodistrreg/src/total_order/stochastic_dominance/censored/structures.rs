use crate::structures::{MAX_OBSERVATIONS, Observation};
use crate::total_order;
use crate::total_order::structures::CovariateStatistic;
use bitree::BITree;

pub type Partition = total_order::structures::Partition<(), ()>;

/// Flag bit marking an observed (uncensored) item in `Estimates::walk_x`.
pub const WALK_OBSERVED_BIT: u32 = 1 << 31;

impl Partition {
    pub fn new(index: usize) -> Self {
        Self {
            index,
            weight: (),
            value: (),
        }
    }
}

/// Combinatorial oracle for "is the interval Kaplan-Meier survival over the covariate
/// interval `[r, s]` exactly 0 once observations `..=data_index` are consumed?", in O(1)
/// per query.
///
/// Derivation: an interval Kaplan-Meier survival is exactly 0 iff the interval's last
/// positive-weight observation (response order, events before censorings at ties) is an
/// EVENT — only then is some factor exactly `1 - w/w`. Clipping preserves this: for any
/// split, the side containing that observation is itself exactly 0 by induction over the
/// span, so the clip interval's lower edge is 0. The condition is purely combinatorial —
/// immune to f32 drift — and monotone in `data_index`: once the interval's last
/// observation is consumed, it holds forever. CDF values that are mathematically exactly
/// 1.0 must come out as exactly 1.0 (downstream consumers like `prediction::mean`'s
/// proper-CDF gate test `== 1.0`), which the estimator cells achieve by consulting this
/// oracle on every write.
///
/// Representation: per covariate `x`, `m[x] = (t[x] << 1) | e[x]`, where `t[x]` is the
/// index into the response-sorted observations of x's last observation and `e[x]` is
/// whether that observation is an event. The `t` values are distinct across covariates
/// (each observation index belongs to exactly one covariate), so a range-max over `m`
/// is achieved by the max-`t` covariate and a single value yields both "which
/// observation is the interval's last" and "is it an event". A covariate that keeps no
/// observations (all censored below every event) retains the initial `m[x] = 0`, which
/// reads as "censored observation at index 0" — even, so it never completes anything,
/// which is exactly right for an empty covariate; the collision with a real censored
/// index-0 observation is harmless for the same reason (and unreachable besides:
/// preprocessing makes the first observation an event).
///
/// Every hot query site grows its interval one covariate at a time (nested-extension
/// order), so the range-max needs no index structure: callers fold [`Self::marker`] into
/// a running max and test it with [`Self::completes_marker`] — O(1) per query, and the
/// oracle stores nothing beyond `m` itself. The linear-scan
/// [`Self::completes_with_all_data`] and [`Self::range_max`] serve the once-per-call
/// finalization path and debug asserts.
#[derive(Debug)]
pub struct CompletionIndex {
    /// `m[x] = (t[x] << 1) | e[x]` as documented on the struct; `m[x] == 0` (no
    /// observation, or a censored observation at index 0) is even and so correctly
    /// never completes.
    m: Vec<u32>,
}

impl CompletionIndex {
    pub fn new(observations: &[Observation<usize, usize, bool, f32>], n_covariate: usize) -> Self {
        // `t` is stored in the upper 31 bits of a u32, which `MAX_OBSERVATIONS` guarantees
        // suffices (the largest index is `len - 1`).
        debug_assert!(observations.len() <= MAX_OBSERVATIONS);

        // Last observation per covariate: index and event bit.
        let mut m = vec![0u32; n_covariate];
        for (index, observation) in observations.iter().enumerate() {
            debug_assert!(observation.weight > 0.0);
            m[observation.x] = ((index as u32) << 1) | u32::from(observation.observed);
        }
        Self { m }
    }

    /// The marker `m[x]` of covariate `x`, for the hot callers' running range-max.
    #[inline(always)]
    pub fn marker(&self, x: usize) -> u32 {
        self.m[x]
    }

    /// Is the interval Kaplan-Meier survival exactly 0 once observations `..=data_index`
    /// are consumed, for an interval whose marker range-max is `marker_max`? True iff
    /// the interval's last observation is an event that has already been consumed.
    /// (`marker_max == 0` ⇒ no observation ⇒ false via the event-bit check.)
    #[inline(always)]
    pub fn completes_marker(marker_max: u32, data_index: usize) -> bool {
        (marker_max & 1) == 1 && ((marker_max >> 1) as usize) <= data_index
    }

    /// Inclusive range-max of `m` over covariates `[r, s]` by linear scan — for the cold
    /// finalization path and debug asserts; hot callers maintain the max incrementally
    /// via [`Self::marker`].
    pub(crate) fn range_max(&self, r: usize, s: usize) -> u32 {
        debug_assert!(r <= s);
        self.m[r..=s].iter().copied().max().unwrap()
    }

    /// [`Self::completes_marker`] with every observation consumed: true iff the
    /// interval's last observation is an event. Linear scan — the only caller walks
    /// disjoint partition ranges once per algorithm call (O(C) total).
    pub fn completes_with_all_data(&self, r: usize, s: usize) -> bool {
        self.range_max(r, s) & 1 == 1
    }
}

/// The Kaplan-Meier running state of one (r, s) cell — only touched when `update_value`
/// actually folds observations, unlike the bounds, which every `propagate_bounds` call
/// rewrites. Kept in a separate array (12 bytes vs 8 for the bounds pair) so each path only
/// streams the bytes it uses.
#[derive(Clone, Debug)]
pub struct SurvivalComputationCold {
    /// The Kaplan-Meier estimator (based on `count` samples). Updated incrementally inside
    /// `update_value`; clipped to the bounds and the result stored as the corresponding entry in
    /// `Estimates::values`.
    pub raw_value: f32,
    /// Sum of weight of included observations.
    pub weight: f32,
    /// Number of samples included in the computation of the raw value. `u32` suffices by
    /// the crate-wide [`MAX_OBSERVATIONS`] bound.
    pub count: u32,
}

/// Lower/upper clip bound of one (r, s) cell — NAN for diagonal singletons. Written as a pair
/// by every `propagate_bounds` call; read as a pair by `update_value`.
pub type Bounds = (f32, f32);

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
    /// Row-major mirror of `values`: `(r, s)` at `compute_row_index((r, s), len)`, so `(r, k)`
    /// for `k` increasing is sequential. The bound-propagation kernel reads a row operand and
    /// a column operand per cell; without the mirror the row side needs a strided gather
    /// (one cache line per element), which measured at ~26% of kernel-bound fits. Every
    /// value write maintains both copies (the second store is sequential and cheap).
    pub values_row: Vec<f32>,
    /// (lower, upper) clip bound of each (r, s). Indexed identically to `values`.
    pub bounds: Vec<Bounds>,
    /// K-M running state of each (r, s). Indexed identically to `values`.
    pub cold: Vec<SurvivalComputationCold>,
    /// Exact-0 completion oracle. Consulted on every cell write so that intervals whose
    /// survival is mathematically exactly 0 are pinned to exactly 0 the moment they
    /// complete, immune to the f32 drift of the running Kaplan-Meier bookkeeping.
    pub completion: CompletionIndex,
    /// Covariate index of each observation (response order, parallel to the context's
    /// `observations`), with [`WALK_OBSERVED_BIT`] set on uncensored items. The hot K-M
    /// walk in `update_value` filters on this pure 4-byte stream instead of dragging
    /// 24-byte `Observation`s through cache; the weight in `walk_w` is only read for
    /// items that pass the filter.
    pub walk_x: Vec<u32>,
    /// Weight of each observation, parallel to `walk_x`.
    pub walk_w: Vec<f32>,
}

impl Estimates {
    pub fn new(
        n: usize,
        start_count: usize,
        completion: CompletionIndex,
        observations: &[Observation<usize, usize, bool, f32>],
    ) -> Self {
        let total = n * (n + 1) / 2;
        let values = vec![1.0; total];
        let values_row = vec![1.0; total];
        let mut bounds = vec![(1.0f32, 1.0f32); total];
        let cold = vec![
            SurvivalComputationCold {
                raw_value: 1.0,
                weight: 0.0,
                count: start_count as u32,
            };
            total
        ];
        // Diagonal entries don't have meaningful bounds (no subdivisions); mark with NAN so the
        // bounds-equality short-circuit in `update_value` would only fire on actually-collapsed
        // intervals.
        for i in 0..n {
            let idx = Self::compute_index((i, i), n);
            bounds[idx] = (f32::NAN, f32::NAN);
        }
        // Covariate indices share their u32 with the flag in the top bit; `MAX_OBSERVATIONS`
        // keeps them below it (the largest covariate index is at most `observations.len() - 1`).
        debug_assert!(observations.len() <= MAX_OBSERVATIONS);
        let walk_x = observations
            .iter()
            .map(|o| o.x as u32 | if o.observed { WALK_OBSERVED_BIT } else { 0 })
            .collect();
        let walk_w = observations.iter().map(|o| o.weight).collect();
        Self {
            len: n,
            values,
            values_row,
            bounds,
            cold,
            completion,
            walk_x,
            walk_w,
        }
    }
    /// Index into the array - `r` is inclusive, `s` is inclusive
    #[inline(always)]
    pub fn compute_index((r, s): (usize, usize), len: usize) -> usize {
        debug_assert!(r <= s);
        debug_assert!(s < len);

        s * (s + 1) / 2 + r
    }
    /// Index into the row-major mirror `values_row` — row `r` (length `len - r`) starts after
    /// rows `0..r`, whose lengths sum to `r * (2 * len + 1 - r) / 2` (the product is always
    /// even, and `r < len` keeps the inner term positive).
    #[inline(always)]
    pub fn compute_row_index((r, s): (usize, usize), len: usize) -> usize {
        debug_assert!(r <= s);
        debug_assert!(s < len);

        r * (2 * len + 1 - r) / 2 + (s - r)
    }
    pub fn len(&self) -> usize {
        self.len
    }

    /// Read the clipped value at (r, s).
    #[inline(always)]
    pub fn value(&self, r: usize, s: usize) -> f32 {
        self.values[Self::compute_index((r, s), self.len)]
    }
    /// Mutable access to the clip bounds at (r, s).
    #[inline(always)]
    pub fn bounds_mut(&mut self, r: usize, s: usize) -> &mut Bounds {
        &mut self.bounds[Self::compute_index((r, s), self.len)]
    }
    /// Write the clipped value at `idx` = (r, s), maintaining the row-major mirror.
    #[inline(always)]
    pub fn set_value(&mut self, idx: usize, row_idx: usize, value: f32) {
        self.values[idx] = value;
        self.values_row[row_idx] = value;
    }
}

impl Estimates {
    /// Bridge from the uncensored-prefix `accelerated_pava` into the generalized-PAVA hot loop.
    /// Both phases run in f32; this just unpacks the BITree's prefix sums into a plain Vec.
    pub fn from_partial_uncensored_solution(
        consumed_weight: BITree<f32>,
        covariate_statistics: &[CovariateStatistic],
        start_count: usize,
        completion: CompletionIndex,
        observations: &[Observation<usize, usize, bool, f32>],
    ) -> Self {
        let len = consumed_weight.len();
        debug_assert_eq!(covariate_statistics.len(), len);

        let mut estimates = Estimates::new(len, start_count, completion, observations);

        let consumed_weight_plain: Vec<f32> = Vec::from(consumed_weight);

        let mut index = 0;
        for s in 0..len {
            let weight_consumed = consumed_weight_plain[s];

            // Filling in the rest of the triangle. The loop runs r descending so the
            // completion queries can fold `marker(r)` into a running range-max over
            // `[r, s]`; the body is order-independent (iteration r only reads
            // `prev_values[r]`/`prev_cold[r]` and writes `curr_values[r]`/`curr_cold[r]`).
            let (prev_values, curr_values) = estimates.values[index - s..].split_at_mut(s);
            let (prev_cold, curr_cold) = estimates.cold[index - s..].split_at_mut(s);
            let mut marker_max = estimates.completion.marker(s);
            for r in (0..s).rev() {
                marker_max = marker_max.max(estimates.completion.marker(r));
                debug_assert_eq!(marker_max, estimates.completion.range_max(r, s));
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
                // Bridge values come from f32 ratio/cumulative-weight arithmetic and a
                // bridged partition can survive untouched to the final row, so completed
                // intervals are pinned to exactly 0 here too. The uncensored prefix is
                // all-observed (`accelerated_pava` stops at the first censored
                // observation), so whenever the interval's last observation lies inside
                // the consumed prefix its event bit holds automatically.
                if CompletionIndex::completes_marker(marker_max, start_count - 1) {
                    curr_cold.raw_value = 0.0;
                }
                *curr_value = curr_cold.raw_value;
            }
            index += s;

            // Diagonal item last; pinned to exactly 0 when complete, like the
            // off-diagonal cells above (`weight_consumed` is a BIT-reconstructed f32 sum
            // and need not equal the covariate's total weight bitwise).
            let diag_raw = if CompletionIndex::completes_marker(
                estimates.completion.marker(s),
                start_count - 1,
            ) {
                0.0
            } else {
                1.0 - weight_consumed / covariate_statistics[s].weight
            };
            let diag_cold = &mut estimates.cold[index];
            diag_cold.raw_value = diag_raw;
            diag_cold.weight = weight_consumed;
            estimates.values[index] = diag_cold.raw_value;
            index += 1;
        }

        // Mirror the freshly built triangle into the row-major copy.
        for r in 0..len {
            for s in r..len {
                let v = estimates.values[Estimates::compute_index((r, s), len)];
                estimates.values_row[Estimates::compute_row_index((r, s), len)] = v;
            }
        }

        estimates
    }
}
