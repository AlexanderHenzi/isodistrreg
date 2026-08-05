use crate::Float;
use crate::routines::binary_search_by_index;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Compute the mean from a (CDF, threshold) pair.
///
/// CDF values are f32 (matching the algorithm's output precision). Thresholds and the
/// returned value live in the caller's response precision `Y`. The integral is accumulated
/// in f64 internally so accumulated round-off stays at f64 precision; we narrow back to `Y`
/// at the return. With ~10⁵ thresholds we'd otherwise lose ~3 f32 ulps to running-sum
/// round-off.
pub fn mean<Y: Float>(
    cdf: impl IntoIterator<Item = f32>,
    thresholds: impl IntoIterator<Item = Y>,
) -> Y {
    let (total, last_cdf_value) = cdf.into_iter().zip(thresholds).fold(
        (0.0f64, 0.0f64),
        |(total, previous), (cdf_value, threshold)| {
            let cdf_value = cdf_value as f64;
            let jump = cdf_value - previous;
            (total + jump * threshold.to_f64().unwrap(), cdf_value)
        },
    );

    // A sub-CDF (censored mass beyond the largest threshold) has an undefined mean. The
    // comparison is exact on purpose: every producer pins values that are mathematically
    // exactly 1 — the uncensored kernels write literal 1.0 rows, `empirical_cdf` pins its
    // final value, and the censored Kaplan-Meier producers pin survival to exactly 0
    // whenever the last positive-weight observation is an event (a purely combinatorial
    // condition, immune to f32 round-off).
    if last_cdf_value < 1.0 {
        Y::nan()
    } else {
        Y::from(total).unwrap()
    }
}

pub trait CovariateInterpolator {
    fn interpolate(&self, response: ResponseCoordinate) -> f32 {
        match response {
            ResponseCoordinate::StrictlyBelowAll => 0.0,
            ResponseCoordinate::AboveOrAtIndex(index) => self.interpolate_index(index),
        }
    }
    fn interpolate_index(&self, index: usize) -> f32;
    fn iter(&self) -> impl ExactSizeIterator<Item = f32> {
        (0..self.len()).map(|i| self.interpolate_index(i))
    }
    fn is_empty(&self) -> bool;
    fn len(&self) -> usize;
}
impl<T: CovariateInterpolator> CovariateInterpolator for &T {
    fn interpolate_index(&self, index: usize) -> f32 {
        (*self).interpolate_index(index)
    }

    fn is_empty(&self) -> bool {
        (*self).is_empty()
    }

    fn len(&self) -> usize {
        (*self).len()
    }
}

pub trait IntoCdfIterator: IntoIterator<Item = f32, IntoIter: ExactSizeIterator> {}
impl<T: IntoIterator<Item = f32, IntoIter: ExactSizeIterator>> IntoCdfIterator for T {}

// Owns the interpolator
pub struct CdfInterpolation<I> {
    interpolator: I,
    index: usize,
    len: usize,
}
impl<I: CovariateInterpolator> CdfInterpolation<I> {
    pub fn new(interpolator: I) -> Self {
        let len = interpolator.len();
        Self {
            interpolator,
            index: 0,
            len,
        }
    }
}
impl<I: CovariateInterpolator> Iterator for CdfInterpolation<I> {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        if self.index < self.len {
            let val = self.interpolator.interpolate_index(self.index);
            self.index += 1;
            Some(val)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.len - self.index;
        (remaining, Some(remaining))
    }
}
impl<I: CovariateInterpolator> ExactSizeIterator for CdfInterpolation<I> {}

impl<const N: usize> CovariateInterpolator for [f32; N] {
    fn interpolate_index(&self, index: usize) -> f32 {
        self[index]
    }

    fn is_empty(&self) -> bool {
        self.as_ref().is_empty()
    }

    fn len(&self) -> usize {
        self.as_ref().len()
    }
}

#[inline]
pub fn search_response<Y: Float>(target: Y, thresholds: &[Y]) -> ResponseCoordinate {
    // An empty threshold grid (e.g. a fit on fully censored data) represents the
    // sub-CDF that is 0 everywhere; binary_search on an empty slice yields Err(0)
    // and thus StrictlyBelowAll, which interpolates to 0.
    debug_assert!(thresholds.windows(2).all(|w| w[0] < w[1]));

    match thresholds.binary_search_by(|t| t.partial_cmp(&target).unwrap()) {
        Err(0) => ResponseCoordinate::StrictlyBelowAll,
        Ok(index) => ResponseCoordinate::AboveOrAtIndex(index),
        Err(index) => ResponseCoordinate::AboveOrAtIndex(index - 1),
    }
}

/// A stateful searcher that, given sorted response values one at a time,
/// yields the corresponding `ResponseCoordinate` against a set of thresholds.
///
/// Thresholds must be non-empty and strictly increasing.
/// Responses must be given in sorted (non-decreasing) order.
pub struct SortedResponseSearcher<'a, Y> {
    thresholds: &'a [Y],
    idx: usize,
}

impl<'a, Y: Float> SortedResponseSearcher<'a, Y> {
    pub fn new(thresholds: &'a [Y]) -> Self {
        debug_assert!(thresholds.windows(2).all(|w| w[0] < w[1]));

        Self { thresholds, idx: 1 }
    }

    /// Feed the next response value (must be >= all previously fed values)
    /// and get back its coordinate.
    pub fn next_response(&mut self, target: Y) -> ResponseCoordinate {
        if self.thresholds.is_empty() || target < self.thresholds[0] {
            ResponseCoordinate::StrictlyBelowAll
        } else {
            while self.idx < self.thresholds.len() && target >= self.thresholds[self.idx] {
                self.idx += 1;
            }

            ResponseCoordinate::AboveOrAtIndex(self.idx - 1)
        }
    }

    /// Wrap a sorted iterator of response values, producing an iterator
    /// of `ResponseCoordinate`s.
    pub fn iter_over<I>(self, responses: I) -> SortedResponseSearch<'a, Y, I::IntoIter>
    where
        I: IntoIterator<Item = Y>,
    {
        SortedResponseSearch {
            searcher: self,
            responses: responses.into_iter(),
        }
    }
}

/// An iterator that feeds sorted response values into a [`SortedResponseSearcher`]
/// and yields the resulting [`ResponseCoordinate`]s.
pub struct SortedResponseSearch<'a, Y, I> {
    searcher: SortedResponseSearcher<'a, Y>,
    responses: I,
}

impl<Y: Float, I> Iterator for SortedResponseSearch<'_, Y, I>
where
    I: Iterator<Item = Y>,
{
    type Item = ResponseCoordinate;

    fn next(&mut self) -> Option<Self::Item> {
        self.responses
            .next()
            .map(|target| self.searcher.next_response(target))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.responses.size_hint()
    }
}

impl<Y: Float, I: ExactSizeIterator<Item = Y>> ExactSizeIterator
    for SortedResponseSearch<'_, Y, I>
{
    fn len(&self) -> usize {
        self.responses.len()
    }
}

/// Targets must be sorted.
pub fn search_responses_sorted<Y: Float, I: IntoIterator<Item = Y>>(
    targets: I,
    thresholds: &[Y],
) -> SortedResponseSearch<'_, Y, I::IntoIter> {
    SortedResponseSearcher::new(thresholds).iter_over(targets)
}

/// Computes quantiles (inverse CDF values) on a discrete threshold grid for a given
/// probability level at one covariate target.
///
/// - Interpolate the CDF across covariates at each threshold index (using `find_coordinates`).
/// - Find an index `i` where the interpolated CDF crosses `probability` using a binary search
///   on `thresholds`.
///
/// Tie‑breaking is controlled via `upper`:
/// - If `upper == false`, returns the smallest threshold `x_i` with `F(x_i) >= probability`
///   (left‑continuous quantile on the grid).
/// - If `upper == true` and `F(x_i) == probability` exactly, returns `x_{i+1}` if it exists;
///   otherwise returns `+∞`.
///
/// If `probability` exceeds all attainable CDF values on the grid (e.g., the last CDF is < 1),
/// `+∞` is returned.
///
/// # Arguments
///
/// - `interpolator`: A `CovariateInterpolator` yielding the f32 CDF values.
/// - `probability`: Probability level in `[0, 1]`. f32 because CDF values are already f32.
/// - `upper`: Tie‑breaking policy when `F(x_i) == probability` (see above).
/// - `thresholds`: Strictly increasing response threshold grid (precision `Y`).
///
/// # Examples
///
/// ```rust
/// use isodistrreg::quantile;
///
/// let thresholds = [10.0, 20.0, 30.0];
///
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.0, false, &thresholds), 10.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.2, false, &thresholds), 10.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.5, false, &thresholds), 20.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.8, false, &thresholds), 20.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 1.0, false, &thresholds), 30.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.0, true, &thresholds), 10.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.2, true, &thresholds), 20.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.5, true, &thresholds), 20.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 0.8, true,  &thresholds), 30.0);
/// assert_eq!(quantile(&[0.2, 0.8, 1.0], 1.0, true, &thresholds), 30.0);
/// ```
///
/// ```rust
/// use isodistrreg::quantile;
///
/// let thresholds = [0.0, 10.0, 20.0, 30.0];
/// assert_eq!(quantile(&[0.0, 0.2, 0.8, 1.0], 0.0, false, &thresholds), 10.0);
/// assert_eq!(quantile(&[0.0, 0.2, 0.8, 1.0], 0.0, true, &thresholds), 10.0);
/// ```
///
/// # Panics
///
/// Panics if `probability` is not in `[0, 1]` (NaN included). The flat-CDF and
/// sub-CDF arms below rely on this range to tell "all mass on the grid" apart from
/// "mass beyond the grid", and a NaN would otherwise fail deep inside the binary
/// search with an unhelpful message.
///
/// ```rust,should_panic
/// isodistrreg::quantile(&[0.5, 1.0], 1.5, false, &[10.0, 20.0]);
/// ```
pub fn quantile<Y: Float, I: CovariateInterpolator>(
    interpolator: &I,
    probability: f32,
    upper: bool,
    thresholds: &[Y],
) -> Y {
    assert!(
        (0.0..=1.0).contains(&probability),
        "probability must be in [0, 1], got {probability}",
    );
    let n = thresholds.len();
    // An empty grid (the empty fit: a sub-CDF that stays at 0) has all its mass beyond
    // every threshold; any quantile of it is +∞.
    if n == 0 {
        return Y::infinity();
    }
    let response_index = binary_search_by_index(0, n, upper, |idx| {
        let coordinate = ResponseCoordinate::AboveOrAtIndex(idx);
        let compare_with = interpolator.interpolate(coordinate);
        compare_with.partial_cmp(&probability).unwrap()
    });
    let index = if upper {
        match response_index {
            // The CDF is flat at exactly `probability` through the last grid index. A
            // mathematically precise upper quantile would always return Y::INFINITY here;
            // we follow the scipy convention and return the supremum of the support
            // instead — the first index attaining the terminal value — but only when
            // that value is 1 (all mass on the grid). For a flat terminal value below 1
            // the missing mass lies beyond the grid, so the support supremum is +∞.
            Ok(index) if index == n - 1 => {
                if probability < 1.0 {
                    None
                } else {
                    // F(n-1) == probability, so an Equal index exists and the
                    // lower-bound search finds the first one — the support supremum.
                    binary_search_by_index(0, n, false, |idx| {
                        let coordinate = ResponseCoordinate::AboveOrAtIndex(idx);
                        let compare_with = interpolator.interpolate(coordinate);
                        compare_with.partial_cmp(&probability).unwrap()
                    })
                    .ok()
                }
            }
            // Value found exactly
            Ok(index) => {
                // By the definition of `binary_search_by_index` this is already the highest
                // index with this cdf value
                Some(index + 1)
            }
            // A sub-CDF, represents Y::INFINITY
            Err(index) if index == n => None,
            // Value not found exactly
            Err(index) => Some(index),
        }
    } else {
        match response_index {
            // A mathematically precise lower quantile would always return Y::NEG_INFINITY
            // here, but we follow the scipy convention and return the infimum of the
            // support
            Err(0) => Some(0),
            Ok(0) if probability == 0.0 => {
                // search for the lower bound of the support
                (0..n).position(|idx| {
                    let coordinate = ResponseCoordinate::AboveOrAtIndex(idx);
                    interpolator.interpolate(coordinate) > 0.0
                })
            }
            // Value found exactly
            Ok(index) => {
                // By the definition of `binary_search_by_index` this is already the lowest
                // index with this cdf value
                Some(index)
            }
            // Value not found exactly
            Err(index) if index < n => {
                // By the definition of `binary_search_by_index` this is already the lowest
                // index with this cdf value
                Some(index)
            }
            // A sub-CDF, represents Y::INFINITY
            Err(_self_len) => None,
        }
    };

    index.map_or(Y::infinity(), |index| thresholds[index])
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ResponseCoordinate {
    StrictlyBelowAll,
    AboveOrAtIndex(usize),
}

#[cfg(test)]
mod test {
    use super::*;

    /// Upper quantile when the CDF is flat at exactly `probability` through the last
    /// grid index: at p = 1 the supremum of the support (first index attaining 1);
    /// below 1 the missing mass lies beyond the grid, so +∞ — and q⁻(p) ≤ q⁺(p).
    #[test]
    fn upper_quantile_flat_terminal_value() {
        // Point mass at 5: support supremum, not the last threshold.
        assert_eq!(quantile(&[1.0f32, 1.0], 1.0, true, &[5.0, 10.0]), 5.0);

        // All-zero sub-CDF: all mass beyond the grid.
        let lower = quantile(&[0.0f32, 0.0], 0.0, false, &[1.0f64, 2.0]);
        let upper = quantile(&[0.0f32, 0.0], 0.0, true, &[1.0f64, 2.0]);
        assert!(upper.is_infinite() && upper > 0.0);
        assert!(lower <= upper);
    }

    /// The empty grid (empty fit: sub-CDF ≡ 0) has all mass beyond every threshold;
    /// every quantile is +∞.
    #[test]
    fn quantile_on_empty_grid_is_infinite() {
        let empty: [f32; 0] = [];
        let q = quantile(&empty, 0.5, false, &[] as &[f64]);
        assert!(q.is_infinite() && q > 0.0);
        let q = quantile(&empty, 0.5, true, &[] as &[f64]);
        assert!(q.is_infinite() && q > 0.0);
    }

    /// Proper CDFs end at exactly 1.0 (every producer pins it); any strictly
    /// smaller last value is a sub-CDF with undefined mean.
    #[test]
    fn mean_requires_exact_unit_mass() {
        let m: f64 = mean([0.5f32, 1.0], [0.0f64, 1.0]);
        assert_eq!(m, 0.5);

        let m: f64 = mean([0.5f32, 1.0 - f32::EPSILON], [0.0f64, 1.0]);
        assert!(m.is_nan());
    }
}
