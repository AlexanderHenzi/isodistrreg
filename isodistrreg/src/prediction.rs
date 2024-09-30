use crate::routines::binary_search_by_index;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub fn mean(cdf: impl IntoIterator<Item = f64>, thresholds: impl IntoIterator<Item = f64>) -> f64 {
    let (total, last_cdf_value) = cdf.into_iter().zip(thresholds).fold(
        (0.0, 0.0),
        |(total, previous), (cdf_value, threshold)| {
            let jump = cdf_value - previous;
            (total + jump * threshold, cdf_value)
        },
    );

    if last_cdf_value < 1.0 {
        f64::NAN
    } else {
        total
    }
}

pub trait CovariateInterpolator {
    fn interpolate(&self, response: ResponseCoordinate) -> f64 {
        match response {
            ResponseCoordinate::StrictlyBelowAll => 0.0,
            ResponseCoordinate::AboveOrAtIndex(index) => self.interpolate_index(index),
        }
    }
    fn interpolate_index(&self, index: usize) -> f64;
    fn iter(&self) -> impl ExactSizeIterator<Item = f64> {
        (0..self.len()).map(|i| self.interpolate_index(i))
    }
    fn is_empty(&self) -> bool;
    fn len(&self) -> usize;
}
impl<T: CovariateInterpolator> CovariateInterpolator for &T {
    fn interpolate_index(&self, index: usize) -> f64 {
        (*self).interpolate_index(index)
    }

    fn is_empty(&self) -> bool {
        (*self).is_empty()
    }

    fn len(&self) -> usize {
        (*self).len()
    }
}

pub trait IntoCdfIterator: IntoIterator<Item = f64, IntoIter: ExactSizeIterator> {}
impl<T: IntoIterator<Item = f64, IntoIter: ExactSizeIterator>> IntoCdfIterator for T {}

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
    type Item = f64;

    fn next(&mut self) -> Option<f64> {
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

impl<const N: usize> CovariateInterpolator for [f64; N] {
    fn interpolate_index(&self, index: usize) -> f64 {
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
pub fn search_response(target: f64, thresholds: &[f64]) -> ResponseCoordinate {
    assert!(!thresholds.is_empty());
    debug_assert!(thresholds.array_windows().all(|&[l, r]| l < r));

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
pub struct SortedResponseSearcher<'a> {
    thresholds: &'a [f64],
    idx: usize,
}

impl<'a> SortedResponseSearcher<'a> {
    pub fn new(thresholds: &'a [f64]) -> Self {
        debug_assert!(thresholds.array_windows().all(|&[l, r]| l < r));

        Self { thresholds, idx: 1 }
    }

    /// Feed the next response value (must be >= all previously fed values)
    /// and get back its coordinate.
    pub fn next_response(&mut self, target: f64) -> ResponseCoordinate {
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
    pub fn iter_over<I>(self, responses: I) -> SortedResponseSearch<'a, I::IntoIter>
    where
        I: IntoIterator<Item = f64>,
    {
        SortedResponseSearch {
            searcher: self,
            responses: responses.into_iter(),
        }
    }
}

/// An iterator that feeds sorted response values into a [`SortedResponseSearcher`]
/// and yields the resulting [`ResponseCoordinate`]s.
pub struct SortedResponseSearch<'a, I> {
    searcher: SortedResponseSearcher<'a>,
    responses: I,
}

impl<I> Iterator for SortedResponseSearch<'_, I>
where
    I: Iterator<Item = f64>,
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

impl<I: ExactSizeIterator<Item = f64>> ExactSizeIterator for SortedResponseSearch<'_, I> {
    fn len(&self) -> usize {
        self.responses.len()
    }
}

/// Targets must be sorted.
pub fn search_responses_sorted<I: IntoIterator<Item = f64>>(
    targets: I,
    thresholds: &[f64],
) -> SortedResponseSearch<'_, I::IntoIter> {
    SortedResponseSearcher::new(thresholds).iter_over(targets)
}

/// Computes quantiles (inverse CDF values) on a discrete threshold grid for given
/// probability levels at one or more covariate targets.
///
/// For each target covariate:
/// - Interpolate the CDF across covariates at each threshold index (using `find_coordinates`).
/// - For each requested probability `q`, find an index `i` where the interpolated CDF
///   crosses `q` using a binary search on `thresholds`.
///
/// Tie‑breaking is controlled via `upper`:
/// - If `upper == false`, returns the smallest threshold `x_i` with `F(x_i) >= q`
///   (left‑continuous quantile on the grid).
/// - If `upper == true` and `F(x_i) == q` exactly, returns `x_{i+1}` if it exists;
///   otherwise returns `+∞`.
///
/// If `q` exceeds all attainable CDF values on the grid (e.g., the last CDF is < 1),
/// `+∞` is returned.
///
/// Results are returned in a flat vector grouped by targets (outer) then quantiles (inner).
///
/// Note: The `quantiles` iterator is collected internally for reuse for each target.
///
/// # Arguments
///
/// - `targets`: Iterator of covariate values at which to compute quantiles.
/// - `quantiles`: Iterator of probability levels in `[0, 1]`.
/// - `upper`: Tie‑breaking policy when `F(x_i) == q`:
///   - `false` → return `x_i` (lower/left rule),
///   - `true`  → return `x_{i+1}` if available, else `+∞`.
/// - `covariates`: Strictly increasing covariate grid of length `C`.
/// - `thresholds`: Strictly increasing response threshold grid of length `N`.
/// - `cdfs`: Flattened CDF matrix of length `C * N` in covariate‑major order.
///
/// # Returns
///
/// A vector of length `targets.count() * quantiles.count()` containing quantile values
/// on the `thresholds` grid for each target and probability (outer: targets; inner: quantiles).
///
/// # Panics
///
/// Panics if `covariates.len() * thresholds.len() != cdfs.len()`.
///
/// # Complexity
///
/// Let `T = number of targets`, `Q = number of quantiles`, `C = covariates.len()`,
/// and `N = thresholds.len()`.
/// - Time: `O(T log C + T * Q * log N)`
/// - Space: `O(T * Q)` for the output (plus `O(Q)` for collected quantiles).
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
pub fn quantile<I: CovariateInterpolator>(
    interpolator: &I,
    probability: f64,
    upper: bool,
    thresholds: &[f64],
) -> f64 {
    let n = thresholds.len();
    let response_index = binary_search_by_index(0, n, upper, |idx| {
        let coordinate = ResponseCoordinate::AboveOrAtIndex(idx);
        let compare_with = interpolator.interpolate(coordinate);
        compare_with.partial_cmp(&probability).unwrap()
    });
    let index = if upper {
        match response_index {
            // A mathematically precise upper quantile would always return f64::INFINITY here,
            // but we follow the scipy convention and return the supremum of the support
            Ok(index) if index == n - 1 => Some(n - 1),
            // Value found exactly
            Ok(index) => {
                // By the definition of `binary_search_by_index` this is already the highest
                // index with this cdf value
                Some(index + 1)
            }
            // A sub-CDF, represents f64::INFINITY
            Err(index) if index == n => None,
            // Value not found exactly
            Err(index) => Some(index),
        }
    } else {
        match response_index {
            // A mathematically precise lower quantile would always return f64::NEG_INFINITY
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
            // A sub-CDF, represents f64::INFINITY
            Err(_self_len) => None,
        }
    };

    index.map_or(f64::INFINITY, |index| thresholds[index])
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ResponseCoordinate {
    StrictlyBelowAll,
    AboveOrAtIndex(usize),
}
