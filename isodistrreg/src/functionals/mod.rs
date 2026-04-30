use crate::structures::Observation;
use std::cmp::Ordering;
#[cfg(feature = "partial-order")]
use std::collections::HashMap;
use std::ops::Range;

mod recursive_clipping;
#[cfg(feature = "partial-order")]
use crate::partial_order::BitSet;
pub use recursive_clipping::ClippingWrapper;

pub trait Functional: Send + Sync {
    type Censoring: Clone + Copy + Default + Send + Sync;
    type Response: Clone + Copy + Default + Send + Sync;
    fn evaluate<Cov>(
        &self,
        obs: impl IntoIterator<Item = Observation<Cov, Self::Response, Self::Censoring>> + Clone,
    ) -> f64;
}

pub trait CauchyMeanValueFunctional {
    type Censoring: Clone + Copy + Default + Send + Sync;
    type Response: Clone + Copy + Default + Send + Sync;

    /// Per-block state held in the algorithm's PAV stack.
    type Block;
    /// Aggregate state shared across blocks (e.g., a clipping triangle). `()` if unused.
    type Shared;

    /// Allocate aggregate state for an n-group problem. Group indices used in
    /// `singleton_block` and `merge_blocks` must lie in `0..n_groups`.
    fn init_shared(&self, n_groups: usize) -> Self::Shared;

    /// Build the block for the singleton range `{group}`.
    fn singleton_block<G, I>(
        &self,
        group: usize,
        get_data: &G,
        shared: &mut Self::Shared,
    ) -> Self::Block
    where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone;

    /// Merge two adjacent blocks, updating `left` in place and consuming `right`.
    /// Pre-condition: `left_range.end == right_range.start`.
    fn merge_blocks<G, I>(
        &self,
        left: &mut Self::Block,
        left_range: Range<usize>,
        right: Self::Block,
        right_range: Range<usize>,
        get_data: &G,
        shared: &mut Self::Shared,
    ) where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone;

    /// Read the current value of a block (used by the algorithm for ordering checks).
    fn block_value(&self, block: &Self::Block) -> f64;

    /// Convenience: evaluate over a contiguous range by sequentially merging singletons.
    /// Used by external callers; the PAV algorithm drives the primitives directly.
    fn evaluate_total_order<G, I>(&self, elements: Range<usize>, get_data: &G) -> f64
    where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    {
        assert!(!elements.is_empty());
        let mut shared = self.init_shared(elements.end);
        let mut acc = self.singleton_block(elements.start, get_data, &mut shared);
        let mut acc_range = elements.start..elements.start + 1;
        for i in (elements.start + 1)..elements.end {
            let next = self.singleton_block(i, get_data, &mut shared);
            self.merge_blocks(
                &mut acc,
                acc_range.clone(),
                next,
                i..i + 1,
                get_data,
                &mut shared,
            );
            acc_range.end = i + 1;
        }
        self.block_value(&acc)
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
        cache: &mut HashMap<BitSet<B>, f64>,
    ) -> f64;
}

/// Marker trait: Can this functional be computed on a single value?
///
/// TODO: Generalize to an integer trait parameter that specifies how many samples are needed?
pub trait SingletonDefinedFunctional: Functional {}

pub struct Average;
impl Average {
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}
impl Default for Average {
    fn default() -> Self {
        Self
    }
}
impl Functional for Average {
    type Censoring = ();
    type Response = f64;
    fn evaluate<Cov>(
        &self,
        obs: impl IntoIterator<Item = Observation<Cov, Self::Response, Self::Censoring>>,
    ) -> f64 {
        let mut count = 0;
        let mut sum = 0.0;
        let mut weight = 0.0;
        obs.into_iter().for_each(|o| {
            count += 1;
            sum += o.weight * o.y;
            weight += o.weight;
        });
        assert!(count >= 1, "need at least one item to compute an average");
        sum / weight
    }
}
impl SingletonDefinedFunctional for Average {}
impl CauchyMeanValueFunctional for Average {
    type Censoring = ();
    type Response = f64;

    /// `(weight, weighted_sum)`. Mean is the quotient.
    type Block = (f64, f64);
    type Shared = ();

    fn init_shared(&self, _n_groups: usize) -> Self::Shared {}

    fn singleton_block<G, I>(
        &self,
        group: usize,
        get_data: &G,
        _shared: &mut Self::Shared,
    ) -> Self::Block
    where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    {
        let (mut weight, mut weighted_sum) = (0.0, 0.0);
        let mut count = 0;
        for o in get_data(group) {
            count += 1;
            weight += o.weight;
            weighted_sum += o.weight * o.y;
        }
        assert!(count >= 1, "need at least one item to compute an average");
        (weight, weighted_sum)
    }

    fn merge_blocks<G, I>(
        &self,
        left: &mut Self::Block,
        _left_range: Range<usize>,
        right: Self::Block,
        _right_range: Range<usize>,
        _get_data: &G,
        _shared: &mut Self::Shared,
    ) where
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    {
        left.0 += right.0;
        left.1 += right.1;
    }

    fn block_value(&self, block: &Self::Block) -> f64 {
        block.1 / block.0
    }

    #[cfg(feature = "partial-order")]
    fn evaluate_partial_order<
        const B: usize,
        G: Fn(usize) -> I + Clone,
        I: Iterator<Item = Observation<(), Self::Response, Self::Censoring>> + Clone,
    >(
        &self,
        elements: &BitSet<B>,
        _successors_inclusive: &[BitSet<B>],
        _predecessors: &[BitSet<B>],
        get_data: &G,
        _cache: &mut HashMap<BitSet<B>, f64>,
    ) -> f64 {
        self.evaluate(elements.iter().flat_map(get_data))
    }
}

pub struct Variance;
impl Variance {
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}
impl Default for Variance {
    fn default() -> Self {
        Self
    }
}
impl Functional for Variance {
    type Censoring = ();
    type Response = f64;
    fn evaluate<Cov>(
        &self,
        obs: impl IntoIterator<Item = Observation<Cov, Self::Response, Self::Censoring>> + Clone,
    ) -> f64 {
        let mut count = 0;
        let (mut inner_product, mut total_weight) = (0.0, 0.0);
        obs.clone().into_iter().for_each(|o| {
            count += 1;
            inner_product += o.weight * o.y;
            total_weight += o.weight;
        });
        if count < 2 {
            // need at least two items to compute a variance
            return f64::NAN;
        }
        let mean = inner_product / total_weight;

        let (mut numerator, mut weight_squared) = (0.0, 0.0);
        obs.into_iter().for_each(|o| {
            numerator += o.weight * (o.y - mean).powi(2);
            weight_squared += o.weight.powi(2);
        });
        let denominator = total_weight - weight_squared / total_weight;

        numerator / denominator
    }
}

pub struct ExpectedShortfall {
    alpha: f64,
}
impl ExpectedShortfall {
    #[must_use]
    pub fn new(alpha: f64) -> Self {
        assert!(alpha > 0.0 && alpha < 1.0, "alpha must be in (0,1)");
        Self { alpha }
    }
}
impl Functional for ExpectedShortfall {
    type Censoring = ();
    type Response = f64;
    fn evaluate<Cov>(
        &self,
        obs: impl IntoIterator<Item = Observation<Cov, Self::Response, Self::Censoring>>,
    ) -> f64 {
        // 1) copy & sort the slice by response ascending
        let mut sorted: Vec<_> = obs.into_iter().collect();
        if sorted.len() < (1.0 / self.alpha).ceil() as usize {
            // need at least this many observations to compute expected shortfall
            return f64::NAN;
        }
        sorted.sort_by(|a, b| {
            // unwrap is safe if you know no NaNs
            a.y.partial_cmp(&b.y).unwrap()
        });

        // 2) compute total weight S and cutoff W* = alpha * S
        let total_w: f64 = sorted.iter().map(|o| o.weight).sum();
        let cutoff: f64 = self.alpha * total_w;

        // 3) walk until we exceed cutoff, accumulating weighted sum
        let mut cum_w = 0.0;
        let mut weighted_sum = 0.0;

        for obs in sorted {
            let w = obs.weight;
            let x = obs.y;

            if cum_w + w < cutoff {
                // take the whole observation
                weighted_sum += w * x;
                cum_w += w;
                continue;
            }
            // else: we only need a partial weight to hit cutoff exactly
            let delta = cutoff - cum_w; // 0 < δ <= w
            weighted_sum += delta * x;
            // we're done
            break;
        }

        // 4) return the weighted‐average of the lower α‐mass
        weighted_sum / cutoff
    }
}

pub trait TotalCmp: Copy + PartialOrd {
    fn total_cmp(&self, other: &Self) -> Ordering;
}

impl TotalCmp for f64 {
    #[inline]
    fn total_cmp(&self, other: &Self) -> Ordering {
        f64::total_cmp(self, other)
    }
}

impl TotalCmp for usize {
    #[inline]
    fn total_cmp(&self, other: &Self) -> Ordering {
        usize::cmp(self, other)
    }
}

#[derive(Debug)]
pub struct KaplanMeier<T> {
    threshold: T,
}
impl<T> KaplanMeier<T> {
    pub fn new(threshold: T) -> Self {
        Self { threshold }
    }
}
impl<T: Default + Sync + Send + TotalCmp> Functional for KaplanMeier<T> {
    type Censoring = bool;
    type Response = T;
    fn evaluate<Cov>(
        &self,
        obs: impl IntoIterator<Item = Observation<Cov, Self::Response, Self::Censoring>>,
    ) -> f64 {
        let mut at_risk = 0.0;
        let mut subset: Vec<_> = obs
            .into_iter()
            .inspect(|o| at_risk += o.weight)
            .filter(|o| o.y <= self.threshold)
            .filter(|o| o.weight > 0.0)
            .map(|o| Observation {
                x: (),
                y: o.y,
                observed: o.observed,
                weight: o.weight,
            })
            .collect();

        subset.sort_unstable_by(|l, r| {
            l.y.total_cmp(&r.y)
                .then(l.observed.cmp(&r.observed).reverse())
        });

        let mut survival = 1.0;
        for o in subset {
            if o.observed {
                survival *= 1.0 - o.weight / at_risk;
            }
            at_risk -= o.weight;
        }

        1.0 - survival
    }
}
impl<T: Default + Sync + Send + TotalCmp> SingletonDefinedFunctional for KaplanMeier<T> {}

#[cfg(test)]
mod test {
    use crate::functionals::{Functional, KaplanMeier};
    use crate::structures::Observation;

    #[test]
    fn kaplan_meier() {
        let km = KaplanMeier::new(3.0);
        assert_eq!(
            km.evaluate(
                [
                    (f64::NAN, 1.0, true),
                    (f64::NAN, 2.0, true),
                    (f64::NAN, 3.0, true),
                ]
                .into_iter()
                .map(Observation::from)
            ),
            1.0,
        );
        assert_eq!(
            km.evaluate(
                [
                    (f64::NAN, 1.0, true),
                    (f64::NAN, 2.0, true),
                    (f64::NAN, 4.0, true),
                ]
                .into_iter()
                .map(Observation::from)
            ),
            2.0 / 3.0,
        );
        assert_eq!(
            km.evaluate(
                [
                    (f64::NAN, 1.0, true),
                    (f64::NAN, 2.0, true),
                    (f64::NAN, 3.0, false),
                ]
                .into_iter()
                .map(Observation::from)
            ),
            2.0 / 3.0,
        );
        assert_eq!(
            km.evaluate(
                [
                    (f64::NAN, 1.0, true),
                    (f64::NAN, 3.0, true),
                    (f64::NAN, 3.0, false),
                ]
                .into_iter()
                .map(Observation::from)
            ),
            2.0 / 3.0,
        );
        assert_eq!(
            km.evaluate(
                [
                    (f64::NAN, 1.0, true),
                    (f64::NAN, 2.0, true),
                    (f64::NAN, 3.0, false),
                    (f64::NAN, 4.0, true),
                ]
                .into_iter()
                .map(Observation::from)
            ),
            0.5,
        );
    }
}
