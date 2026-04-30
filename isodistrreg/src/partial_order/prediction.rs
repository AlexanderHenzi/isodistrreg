use crate::partial_order::routines::find_neighbors;
use crate::partial_order::{OrderingInfo, PredictionWorkspace};
use crate::prediction::{CdfInterpolation, CovariateInterpolator};
use std::mem;

pub struct Interpolation<'a> {
    lower_neighbors: Vec<usize>,
    upper_neighbors: Vec<usize>,
    cdfs: &'a [f64],
    global_cdf: &'a [f64],
}

impl<'a> Interpolation<'a> {
    pub fn new(
        target: &[f64],
        covariates: &[f64],
        ordering_info: &OrderingInfo,
        increasing: bool,
        (cdfs, global_cdf): (&'a [f64], &'a [f64]),
        workspace: &mut PredictionWorkspace,
    ) -> Self {
        let mut lower_neighbors: Vec<_> =
            find_neighbors(target, covariates, ordering_info, workspace, true).collect();
        let mut upper_neighbors: Vec<_> =
            find_neighbors(target, covariates, ordering_info, workspace, false).collect();

        if !increasing {
            mem::swap(&mut lower_neighbors, &mut upper_neighbors);
        }

        Self {
            lower_neighbors,
            upper_neighbors,
            cdfs,
            global_cdf,
        }
    }

    fn get_cdf_value(&self, i: usize, j: usize) -> f64 {
        let n_threshold = self.global_cdf.len();
        self.cdfs[i * n_threshold + j]
    }
}

impl CovariateInterpolator for Interpolation<'_> {
    fn interpolate_index(&self, index: usize) -> f64 {
        let upper_bound = self
            .lower_neighbors
            .iter()
            .map(|&i| self.get_cdf_value(i, index))
            .reduce(f64::min);
        let lower_bound = self
            .upper_neighbors
            .iter()
            .map(|&i| self.get_cdf_value(i, index))
            .reduce(f64::max);
        match (lower_bound, upper_bound) {
            (Some(lower), Some(upper)) => lower.midpoint(upper),
            (Some(lower), None) => lower,
            (None, Some(upper)) => upper,
            (None, None) => self.global_cdf[index],
        }
    }

    fn is_empty(&self) -> bool {
        self.global_cdf.is_empty()
    }

    fn len(&self) -> usize {
        self.global_cdf.len()
    }
}

impl IntoIterator for Interpolation<'_> {
    type Item = f64;
    type IntoIter = CdfInterpolation<Self>;

    fn into_iter(self) -> Self::IntoIter {
        CdfInterpolation::new(self)
    }
}
