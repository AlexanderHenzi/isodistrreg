use crate::prediction::CdfInterpolation;
use crate::prediction::CovariateInterpolator;

/// Interpolates at a given covariate value.
pub enum Interpolation<'a> {
    /// Covariate lands exactly on an index: read directly from one CDF.
    Exact { cdf: &'a [f64] },
    /// Covariate falls between two indices: linearly blend two CDFs.
    Between {
        left_cdf: &'a [f64],
        right_cdf: &'a [f64],
        share: f64,
    },
}

impl<'a> Interpolation<'a> {
    #[must_use]
    pub fn new(target: f64, covariates: &[f64], cdfs: (&'a [f64], usize)) -> Self {
        Self::from_coordinate(search_covariate(target, covariates), cdfs)
    }
    #[must_use]
    pub fn from_coordinate(target: CovariateCoordinate, cdfs: (&'a [f64], usize)) -> Self {
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
    fn get_cdf(i: usize, (cdfs, n_threshold): (&[f64], usize)) -> &[f64] {
        debug_assert!(!cdfs.is_empty());
        debug_assert_ne!(n_threshold, 0);
        debug_assert_eq!(cdfs.len() % n_threshold, 0);

        &cdfs[i * n_threshold..(i + 1) * n_threshold]
    }
}

impl CovariateInterpolator for Interpolation<'_> {
    fn interpolate_index(&self, index: usize) -> f64 {
        match self {
            Interpolation::Exact { cdf } => cdf[index],
            Interpolation::Between {
                left_cdf,
                right_cdf,
                share,
            } => {
                let left = left_cdf[index];
                let right = right_cdf[index];
                (1.0 - share) * left + share * right
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
    type Item = f64;
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
/// ascending) finite `f64` values.
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
///   `0.0 < share < 1.0`. This makes `target` recoverable as
///   `(1.0 - share) * slice[i] + share * slice[i + 1]`.
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
pub fn search_covariate(target: f64, slice: &[f64]) -> CovariateCoordinate {
    debug_assert!(!slice.is_empty());
    debug_assert!(slice.array_windows().all(|&[l, r]| l < r));

    match slice.binary_search_by(|c| c.partial_cmp(&target).unwrap()) {
        Ok(index) => CovariateCoordinate::Exact(index),
        Err(index) => match index {
            0 => CovariateCoordinate::Exact(0),
            in_range if in_range < slice.len() => {
                let range = slice[index] - slice[index - 1];
                let share = (target - slice[index - 1]) / range;
                CovariateCoordinate::Between {
                    left: index - 1,
                    share,
                }
            }
            after_last => {
                debug_assert_eq!(after_last, slice.len());
                CovariateCoordinate::Exact(slice.len() - 1)
            }
        },
    }
}

pub struct GridPredictorState<'a> {
    pub search: CovariateSearch<'a>,
    pub interpolation: Interpolation<'a>,
    pub cdfs: (&'a [f64], usize),
}
impl GridPredictorState<'_> {
    pub fn update(&mut self, query: f64) {
        let coordinate = self.search.advance(query);
        self.interpolation = Interpolation::from_coordinate(coordinate, self.cdfs);
    }
}

pub struct CovariateSearch<'a> {
    references: &'a [f64],
    idx: usize,
}

impl<'a> CovariateSearch<'a> {
    #[must_use]
    pub fn new(references: &'a [f64]) -> Self {
        assert!(!references.is_empty());
        Self { references, idx: 0 }
    }

    /// Advance the search with the next query value (must be provided in-order).
    pub fn advance(&mut self, q: f64) -> CovariateCoordinate {
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
                let left = self.idx - 1;
                let right = self.idx;
                let range = self.references[right] - self.references[left];
                let share = (q - self.references[left]) / range;
                CovariateCoordinate::Between { left, share }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CovariateCoordinate {
    /// At this exact index
    Exact(usize),
    /// Between left and left+1, 0.0 < share < 1.0
    Between { left: usize, share: f64 },
}

#[cfg(test)]
mod test {
    use crate::total_order::prediction::{CovariateCoordinate, search_covariate};

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
