use crate::Float;
use crate::prediction::CdfInterpolation;
use crate::prediction::CovariateInterpolator;

/// Interpolates at a given covariate value.
pub enum Interpolation<'a> {
    /// Covariate lands exactly on an index: read directly from one CDF.
    Exact { cdf: &'a [f32] },
    /// Covariate falls between two indices: linearly blend two CDFs.
    Between {
        left_cdf: &'a [f32],
        right_cdf: &'a [f32],
        share: f32,
    },
}

impl<'a> Interpolation<'a> {
    #[must_use]
    pub fn new<X: Float>(target: X, covariates: &[X], cdfs: (&'a [f32], usize)) -> Self {
        Self::from_coordinate(search_covariate(target, covariates), cdfs)
    }
    #[must_use]
    pub fn from_coordinate(target: CovariateCoordinate, cdfs: (&'a [f32], usize)) -> Self {
        if cdfs.0.is_empty() {
            // The empty fit: a sub-CDF that stays at 0.0, represented by an empty CDF
            // slice (every response then resolves to `StrictlyBelowAll` -> 0.0).
            return Self::Exact { cdf: cdfs.0 };
        }
        match target {
            CovariateCoordinate::Exact(index) => Self::Exact {
                cdf: Self::get_cdf(index, cdfs),
            },
            CovariateCoordinate::Between { left, share } => Self::Between {
                left_cdf: Self::get_cdf(left, cdfs),
                right_cdf: Self::get_cdf(left + 1, cdfs),
                share,
            },
        }
    }
    fn get_cdf(i: usize, (cdfs, n_threshold): (&[f32], usize)) -> &[f32] {
        debug_assert!(!cdfs.is_empty());
        debug_assert_ne!(n_threshold, 0);
        debug_assert_eq!(cdfs.len() % n_threshold, 0);

        &cdfs[i * n_threshold..(i + 1) * n_threshold]
    }
}

impl CovariateInterpolator for Interpolation<'_> {
    fn interpolate_index(&self, index: usize) -> f32 {
        match self {
            Interpolation::Exact { cdf } => cdf[index],
            Interpolation::Between {
                left_cdf,
                right_cdf,
                share,
            } => {
                let left = left_cdf[index];
                let right = right_cdf[index];
                // Numerically-stable monotone lerp. `(1 - share) * left + share * right`
                // does not reproduce the common value when `left == right`: for two f32
                // endpoints that are both exactly `p` (e.g. a CDF pinned at a quantile
                // level by Kaplan-Meier), it can round to `p - 1 ulp` for one share and
                // exactly `p` for another. Because the CDF is stochastically ordered in
                // the covariate (`left >= right`), that sub-ulp wobble flips a `F >= p`
                // quantile test between adjacent covariates and breaks monotonicity of the
                // quantile in x. This form has `right - left == 0` exactly when the
                // endpoints agree, so it returns `left` for every share, and it is monotone
                // in `share` for a fixed-sign delta.
                left + share * (right - left)
            }
        }
    }

    fn is_empty(&self) -> bool {
        match self {
            Interpolation::Exact { cdf } => cdf.is_empty(),
            Interpolation::Between {
                left_cdf,
                right_cdf,
                ..
            } => {
                debug_assert_eq!(left_cdf.len(), right_cdf.len());
                left_cdf.is_empty()
            }
        }
    }

    fn len(&self) -> usize {
        match self {
            Interpolation::Exact { cdf } => cdf.len(),
            Interpolation::Between {
                left_cdf,
                right_cdf,
                ..
            } => {
                debug_assert_eq!(left_cdf.len(), right_cdf.len());
                left_cdf.len()
            }
        }
    }
}

impl IntoIterator for Interpolation<'_> {
    type Item = f32;
    type IntoIter = CdfInterpolation<Self>;

    fn into_iter(self) -> Self::IntoIter {
        CdfInterpolation::new(self)
    }
}

/// Locate a target within a sorted slice and, if not an exact match, return the left bounding index
/// together with a linear interpolation weight.
///
/// This helper wraps a binary search over `slice` and interprets the insertion point as a
/// coordinate for interpolation. It is intended for slices of monotonically increasing (strictly
/// ascending) finite floating-point values.
///
/// # Returns
///
/// - `Ok(i)` if `target` is exactly equal to `slice[i]`.
/// - `Ok(0)` if `target` is less than the first element (`target <= slice[0]`).
/// - `Ok(len - 1)` if `target` is greater than the last element
///   (`target >= slice[len - 1]`).
/// - `Err((i, share))` if `target` lies strictly between `slice[i]` and
///   `slice[i + 1]`, where
///   `share = (target - slice[i]) / (slice[i + 1] - slice[i])` and thus
///   `0.0 < share < 1.0`. The `share` is computed in `X` arithmetic to avoid catastrophic
///   cancellation when `X = f64` (two f64-distinct neighbors can collapse to a single f32),
///   then narrowed to f32 for the CDF blend. If the `X`-precision share lands outside the open
///   interval `(0, 1)` (or the endpoints would otherwise be degenerate), the function falls
///   back to an exact endpoint instead.
///
/// # Requirements
///
/// - `slice` must be non-empty.
/// - Values in `slice` must be strictly increasing.
/// - Inputs should be finite (no `NaN`).
///
/// # Complexity
/// - `O(log n)` comparisons, where `n = slice.len()`.
///
/// # Examples
///
/// See the module tests.
#[inline]
#[must_use]
pub fn search_covariate<X: Float>(target: X, slice: &[X]) -> CovariateCoordinate {
    debug_assert!(!slice.is_empty());
    debug_assert!(slice.windows(2).all(|w| w[0] < w[1]));

    match slice.binary_search_by(|c| c.partial_cmp(&target).unwrap()) {
        Ok(index) => CovariateCoordinate::Exact(index),
        Err(index) => match index {
            0 => CovariateCoordinate::Exact(0),
            in_range if in_range < slice.len() => between_or_snap(target, slice, index - 1),
            after_last => {
                debug_assert_eq!(after_last, slice.len());
                CovariateCoordinate::Exact(slice.len() - 1)
            }
        },
    }
}

/// Place `target` between `slice[left]` and `slice[left + 1]`, returning a `Between` with the
/// linear interpolation share or snapping to an exact endpoint when the share is degenerate.
///
/// The share is computed in `X` precision and narrowed to f32 only at the end. This avoids
/// catastrophic cancellation that would occur if we narrowed `target` / `lo` / `hi` to f32
/// first — at large magnitudes two f64-distinct neighbors can round to the same f32, leaving
/// `hi - lo == 0` and producing ±∞ or NaN. The share is used downstream to blend f32 CDF
/// endpoints, so f32 storage of the share itself is sufficient.
///
/// Classification (happy path first):
/// - `0 < share < 1`: blend → `Between { left, share }`.
/// - `share >= 1` (incl. `+∞`): snap to `Exact(left + 1)`.
/// - `share <= 0`, `NaN`, or `-∞`: snap to `Exact(left)`. NaN lands here because every
///   comparison with NaN is `false`, so it fails both `share > 0` (no happy path) and
///   `share >= 1` (no right snap), and falls through to the final arm.
///
/// Callers must ensure `left + 1 < slice.len()`.
#[inline]
fn between_or_snap<X: Float>(target: X, slice: &[X], left: usize) -> CovariateCoordinate {
    let lo = slice[left];
    let hi = slice[left + 1];
    let share_x = (target - lo) / (hi - lo);

    if share_x > X::zero() && share_x < X::one() {
        CovariateCoordinate::Between {
            left,
            share: share_x.to_f32().unwrap(),
        }
    } else if share_x >= X::one() {
        CovariateCoordinate::Exact(left + 1)
    } else {
        CovariateCoordinate::Exact(left)
    }
}

pub struct GridPredictorState<'a, X> {
    pub search: CovariateSearch<'a, X>,
    pub interpolation: Interpolation<'a>,
    pub cdfs: (&'a [f32], usize),
}
impl<X: Float> GridPredictorState<'_, X> {
    pub fn update(&mut self, query: X) {
        let coordinate = self.search.advance(query);
        self.interpolation = Interpolation::from_coordinate(coordinate, self.cdfs);
    }
}

pub struct CovariateSearch<'a, X> {
    references: &'a [X],
    idx: usize,
}

impl<'a, X: Float> CovariateSearch<'a, X> {
    #[must_use]
    pub fn new(references: &'a [X]) -> Self {
        Self { references, idx: 0 }
    }

    /// Advance the search with the next query value (must be provided in-order).
    pub fn advance(&mut self, q: X) -> CovariateCoordinate {
        if self.references.is_empty() {
            // The empty fit has no covariate grid; the coordinate is irrelevant because
            // its (empty) CDF interpolates to 0.0 for every response.
            return CovariateCoordinate::Exact(0);
        }
        let last = self.references.len() - 1;

        if q <= self.references[0] {
            CovariateCoordinate::Exact(0)
        } else if q >= self.references[last] {
            CovariateCoordinate::Exact(last)
        } else {
            while self.idx < self.references.len() && self.references[self.idx] < q {
                self.idx += 1;
            }

            if q == self.references[self.idx] {
                CovariateCoordinate::Exact(self.idx)
            } else {
                between_or_snap(q, self.references, self.idx - 1)
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CovariateCoordinate {
    /// At this exact index
    Exact(usize),
    /// Between left and left+1, 0.0 < share < 1.0
    Between { left: usize, share: f32 },
}

#[cfg(test)]
mod test {
    use crate::prediction::{CovariateInterpolator, quantile};
    use crate::total_order::prediction::{CovariateCoordinate, Interpolation, search_covariate};

    /// The covariate blend must reproduce a shared endpoint value EXACTLY, for every
    /// share. `(1 - share) * v + share * v` does not: `0.9` is not representable in
    /// f32, and that form rounds it to `0.9 - 1 ulp` for some shares while returning
    /// exactly `0.9` for others. Downstream that sub-ulp wobble flips a `F >= 0.9`
    /// quantile test between neighbouring covariates and breaks monotonicity of the
    /// quantile in x (regression test for the covariate interpolation).
    #[test]
    fn blend_reproduces_equal_endpoints_across_shares() {
        // Values deliberately not representable in binary (0.1, 0.9, ...) so a naive
        // blend drifts off them.
        let cdf: [f32; 4] = [0.1, 0.5, 0.9, 1.0];
        for k in 0..=1000 {
            let share = k as f32 / 1000.0;
            let interpolation = Interpolation::Between {
                left_cdf: &cdf,
                right_cdf: &cdf,
                share,
            };
            for index in 0..cdf.len() {
                assert_eq!(
                    interpolation.interpolate_index(index),
                    cdf[index],
                    "share {share} at index {index} did not reproduce the common value",
                );
            }
        }
    }

    /// For stochastically ordered endpoints (`left >= right` pointwise, i.e. a smaller
    /// covariate has the larger CDF), the interpolated CDF must be non-increasing in
    /// `share` at every response index — sweeping the covariate up may only shift mass
    /// upward. A single sub-ulp increase here is what lets a quantile move the wrong way.
    #[test]
    fn blend_pointwise_non_increasing_in_share() {
        let left: [f32; 5] = [0.3, 0.9, 0.9, 0.95, 1.0];
        let right: [f32; 5] = [0.1, 0.5, 0.9, 0.90, 1.0];
        assert!(left.iter().zip(&right).all(|(l, r)| l >= *r));

        let mut previous = [f32::INFINITY; 5];
        for k in 0..=1000 {
            let share = k as f32 / 1000.0;
            let interpolation = Interpolation::Between {
                left_cdf: &left,
                right_cdf: &right,
                share,
            };
            for index in 0..left.len() {
                let value = interpolation.interpolate_index(index);
                assert!(
                    value <= previous[index],
                    "share {share} at index {index}: {value} > {} (increased with share)",
                    previous[index],
                );
                previous[index] = value;
            }
        }
    }

    /// End goal of the two invariants above: the lower quantile is non-decreasing as the
    /// covariate rises. The endpoints share the exact value `0.9` at threshold index 2 —
    /// the adversarial case that motivated the fix — so the pre-fix blend would round it
    /// below `0.9` for some shares and bounce the quantile between `3.0` and `4.0`.
    #[test]
    fn quantile_non_decreasing_along_covariate_blend() {
        let left: [f32; 5] = [0.3, 0.9, 0.9, 0.95, 1.0];
        let right: [f32; 5] = [0.1, 0.5, 0.9, 0.90, 1.0];
        let thresholds = [1.0f64, 2.0, 3.0, 4.0, 5.0];

        let mut previous = f64::NEG_INFINITY;
        for k in 0..=1000 {
            let share = k as f32 / 1000.0;
            let interpolation = Interpolation::Between {
                left_cdf: &left,
                right_cdf: &right,
                share,
            };
            let q = quantile(&interpolation, 0.9, false, &thresholds);
            assert!(
                q >= previous,
                "share {share}: quantile {q} dropped below previous {previous}",
            );
            previous = q;
        }
        // Sanity: the quantile actually moved (endpoints are not identical), so the test
        // exercises a real transition rather than a constant.
        assert!(previous > 2.0);
    }

    #[test]
    fn exact_match() {
        assert_eq!(
            search_covariate(3.0, &[1.0, 3.0, 5.0]),
            CovariateCoordinate::Exact(1)
        );
        assert_eq!(
            search_covariate(0.0, &[-1.0, -0.0, 2.0]),
            CovariateCoordinate::Exact(1)
        );
    }

    #[test]
    fn between() {
        let xs = [0.0, 10.0, 20.0];
        let target = 15.0;
        let result = search_covariate(target, &xs);
        if let CovariateCoordinate::Between { left, share } = result {
            assert_eq!(left, 1);
            assert_eq!(share, 0.5);
            let share = f64::from(share);
            let reconstructed = (1.0 - share) * xs[left] + share * xs[left + 1];
            assert!((reconstructed - target).abs() < 1e-12);
        } else {
            assert!(false);
        }
    }

    #[test]
    fn extremes() {
        assert_eq!(
            search_covariate(-3.0, &[10.0, 20.0, 30.0]),
            CovariateCoordinate::Exact(0)
        );
        assert_eq!(
            search_covariate(100.0, &[10.0, 20.0, 30.0]),
            CovariateCoordinate::Exact(2)
        );
    }
}
