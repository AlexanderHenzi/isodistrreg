use crate::functionals::{CauchyMeanValueFunctional, Functional};
#[cfg(feature = "partial-order")]
use crate::partial_order::BitSet;
use crate::structures::Observation;
use num_traits::Float as NumFloat;
#[cfg(feature = "partial-order")]
use std::collections::HashMap;
use std::ops::Range;

/// Packed-triangle index for cell `(r, s)` with `r <= s`.
#[inline]
fn index(r: usize, s: usize) -> usize {
    s * (s + 1) / 2 + r
}

pub struct ClippingWrapper<F: Functional>(F);
impl<F: Functional> ClippingWrapper<F> {
    pub fn new(functional: F) -> Self {
        Self(functional)
    }
}

impl<F: Functional> CauchyMeanValueFunctional for ClippingWrapper<F> {
    type Censoring = F::Censoring;
    type Response = F::Response;
    type Value = F::Value;

    /// The clipped value for this block's range.
    type Block = F::Value;
    /// Packed lower triangle indexed by absolute group indices: cell `(r, s)` (with `r <= s`)
    /// is at `s * (s + 1) / 2 + r`. Cells live for the duration of one algorithm run.
    type Shared = Vec<F::Value>;

    fn init_shared(&self, n_groups: usize) -> Self::Shared {
        vec![F::Value::nan(); n_groups * (n_groups + 1) / 2]
    }

    fn singleton_block<G, I>(
        &self,
        group: usize,
        get_data: &G,
        shared: &mut Self::Shared,
    ) -> Self::Block
    where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    {
        let value = self.0.evaluate(get_data(group));
        shared[index(group, group)] = value;
        value
    }

    fn merge_blocks<G, I>(
        &self,
        left: &mut Self::Block,
        left_range: Range<usize>,
        _right: Self::Block,
        right_range: Range<usize>,
        get_data: &G,
        shared: &mut Self::Shared,
    ) where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    {
        debug_assert_eq!(left_range.end, right_range.start);

        // Fill the rectangle `(r, s)` with `r ∈ [left.start, left.end]` and
        // `s ∈ [left.end + 1, right.end]`, in increasing length so every dependency lies either in
        // a previously-built sub-triangle (left/right) or in a rectangle cell already filled in an
        // earlier iteration of this loop.
        for length in 2..=((right_range.end - 1) - left_range.start + 1) {
            for r in left_range.start..=left_range.end - 1 {
                let s = r + length - 1;
                if s < (left_range.end - 1) + 1 || s > (right_range.end - 1) {
                    continue;
                }

                let raw = self.0.evaluate((r..=s).flat_map(get_data));

                let (mut lower, mut upper) = (F::Value::neg_infinity(), F::Value::infinity());
                for k in r..s {
                    let l = shared[index(r, k)];
                    let r = shared[index(k + 1, s)];
                    // A split with an undefined side imposes no clipping
                    // constraint, and depending on the functional, the number
                    // of elements needed to compute a value can be > 1, such
                    // as with the sample variance. NaN-ignoring min/max would
                    // instead pin the value to the defined side.
                    if l.is_nan() || r.is_nan() {
                        continue;
                    }
                    lower = lower.max(l.min(r));
                    upper = upper.min(l.max(r));
                }

                shared[index(r, s)] = raw.max(lower).min(upper);
            }
        }

        *left = shared[index(left_range.start, right_range.end - 1)];
    }

    fn block_value(&self, block: &Self::Block) -> Self::Value {
        *block
    }

    #[cfg(feature = "partial-order")]
    fn evaluate_partial_order<
        const B: usize,
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    >(
        &self,
        elements: &BitSet<B>,
        successors_inclusive: &[BitSet<B>],
        predecessors: &[BitSet<B>],
        get_data: &G,
        cache: &mut HashMap<BitSet<B>, F::Value>,
    ) -> F::Value {
        if let Some(result) = cache.get(elements) {
            return *result;
        }

        let unclipped_value = self.0.evaluate(elements.iter().flat_map(get_data));

        let (mut lower_bound, mut upper_bound) = (F::Value::neg_infinity(), F::Value::infinity());
        dfs(
            BitSet::new(),
            &BitSet::new(),
            elements,
            &mut lower_bound,
            &mut upper_bound,
            successors_inclusive,
            predecessors,
            get_data,
            self,
            cache,
        );

        let value = unclipped_value.max(lower_bound).min(upper_bound);

        cache.insert(elements.clone(), value);

        value
    }
}

#[cfg(feature = "partial-order")]
#[allow(clippy::too_many_arguments)]
fn dfs<
    const B: usize,
    F: Functional,
    G: Fn(usize) -> I + Clone,
    I: Iterator<Item = Observation<(), F::Response, F::Censoring>> + Clone,
>(
    in_set: BitSet<B>,
    forbidden: &BitSet<B>,
    universe: &BitSet<B>,
    lower_bound: &mut F::Value,
    upper_bound: &mut F::Value,
    successors_inclusive: &[BitSet<B>],
    predecessors: &[BitSet<B>],
    get_data: &G,
    recursive_functional: &ClippingWrapper<F>,
    cache: &mut HashMap<BitSet<B>, F::Value>,
) {
    let decided = in_set.union(forbidden);
    let undecided = universe.difference(&decided);

    // Available v: undecided and all predecessors (incl) of v are already in in_set
    // (equivalently, strict preds ⊆ in_set since v itself is not in_set yet).
    match undecided.min() {
        None => {
            // We found a lower set that we should check

            // No bounds for these extremes
            let lower_set = in_set;
            if lower_set.is_empty() || lower_set.eq(universe) {
                return;
            }
            let upper_set = universe.difference(&lower_set);

            let lower_set_value = recursive_functional.evaluate_partial_order(
                &lower_set,
                successors_inclusive,
                predecessors,
                get_data,
                cache,
            );
            if lower_set_value.is_nan() {
                return;
            }

            let upper_set_value = recursive_functional.evaluate_partial_order(
                &upper_set,
                successors_inclusive,
                predecessors,
                get_data,
                cache,
            );
            if upper_set_value.is_nan() {
                return;
            }

            let lowest = lower_set_value.min(upper_set_value);
            let highest = lower_set_value.max(upper_set_value);

            *lower_bound = lower_bound.max(lowest);
            *upper_bound = upper_bound.min(highest);
        }
        Some(minimal_element) => {
            // include v
            dfs(
                in_set.with_singleton(minimal_element),
                forbidden,
                universe,
                lower_bound,
                upper_bound,
                successors_inclusive,
                predecessors,
                get_data,
                recursive_functional,
                cache,
            );

            // exclude v => forbid v and everything above it (keeps "lower set" property)
            dfs(
                in_set,
                &forbidden.union(&successors_inclusive[minimal_element]),
                universe,
                lower_bound,
                upper_bound,
                successors_inclusive,
                predecessors,
                get_data,
                recursive_functional,
                cache,
            );
        }
    }
}
