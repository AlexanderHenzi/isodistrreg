use crate::IntoCdfIterator;
use crate::error::Error;
#[cfg(feature = "partial-order")]
use crate::partial_order;
#[cfg(feature = "partial-order")]
use crate::partial_order::PredictionWorkspace;
use crate::prediction::{
    CdfInterpolation, CovariateInterpolator, ResponseCoordinate, search_responses_sorted,
};
use crate::preprocessing::validate;
use crate::routines::{lexicographic_order, transpose};
use crate::structures::StochasticOrder;
use crate::total_order;
use crate::total_order::{
    CovariateSearch, GridPredictorState, Interpolation as SingleInterpolation,
};
use crate::{IsotonicDistributionalRegressionFit, ProgressTracker};
use itertools::Either;
use rand::distr::Distribution;
use rand::distr::Uniform;
use rand::seq::index;
#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelIterator, ParallelIterator};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::sync::Mutex;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Fit<T: IsotonicDistributionalRegressionFit> {
    pub fits: Vec<T>,
    /// All covariates that were available to this problem, deduplicated.
    pub covariates: Vec<f64>,
    /// All uncensored responses that were available to this problem, deduplicated.
    pub thresholds: Vec<f64>,
    /// For each fit, map the global threshold to the fit of the threshold.
    pub threshold_map: Vec<ResponseCoordinate>,
    pub covariate_groups: T::CovariateOrder,
}

impl<T: IsotonicDistributionalRegressionFit> Fit<T> {
    /// Creates a [`Fit`] directly from pre-computed components.
    pub fn from_parts(
        fit: T,
        covariates: Vec<f64>,
        thresholds: Vec<f64>,
        covariate_groups: T::CovariateOrder,
    ) -> Self
    where
        T::CovariateOrder: Sized,
    {
        let threshold_map = (0..thresholds.len())
            .map(ResponseCoordinate::AboveOrAtIndex)
            .collect();

        Self {
            fits: vec![fit],
            covariates,
            thresholds,
            threshold_map,
            covariate_groups,
        }
    }
}

pub struct Config {
    n_subsamples: usize,
    subsample_size: usize,
    replace: bool,
    /// Random seed for reproducibility. Only relevant when (su)bagging is active
    /// (i.e. when more than one subsample is drawn).
    seed: Option<u64>,
    /// Number of worker threads used to fit the individual subsamples in parallel.
    /// Only relevant when (su)bagging is active and the `parallel` feature is
    /// enabled. A value of `1` runs the subfits serially.
    #[cfg_attr(not(feature = "parallel"), allow(dead_code))]
    n_jobs: usize,
}
impl Config {
    #[must_use]
    pub fn new(
        n_subsamples: usize,
        subsample_size: usize,
        replace: bool,
        seed: Option<u64>,
        n_jobs: usize,
    ) -> Self {
        debug_assert!(n_jobs >= 1, "n_jobs must be at least 1");
        Self {
            n_subsamples,
            subsample_size,
            replace,
            seed,
            n_jobs,
        }
    }

    #[must_use]
    pub fn disable(subsample_size: usize) -> Self {
        Self::new(1, subsample_size, false, None, 1)
    }

    pub fn parse(
        subsamples: Option<usize>,
        subsample_size: Option<Either<usize, f64>>,
        replace: bool,
        n: usize,
        seed: Option<u64>,
        n_jobs: usize,
    ) -> Result<Config, Error> {
        debug_assert_ne!(n, 0);

        use crate::Error::SubaggingParameterInconsistency as BadParam;
        use Either::{Left, Right};

        // Validate individual values
        if matches!(subsamples, Some(0)) {
            return Err(BadParam("subsample count should be >= 1"));
        }
        if matches!(subsample_size, Some(Left(c)) if c == 0 || c > n) {
            return Err(BadParam(
                "subsample size should be at least 1 and at most n if provided as integer",
            ));
        }
        if matches!(subsample_size, Some(Right(s)) if !(0.0 < s && s <= 1.0)) {
            return Err(BadParam(
                "subsample size should be in (0.0, 1.0] if provided as float",
            ));
        }
        if n_jobs == 0 {
            return Err(BadParam("n_jobs must be at least 1"));
        }

        // Simplify by recognizing default values
        let subsample_size = match subsample_size {
            Some(Left(count)) if count == n => None,
            Some(Right(1.0)) => None,
            other => other,
        };

        match (subsamples, subsample_size) {
            (None, Some(_)) => Err(BadParam(
                "specify the number of subsamples when specifying a subsample size",
            )),
            (None, None) if replace => Err(BadParam(
                "can't sample with replacement if subsample_size not specified",
            )),
            // No (su)bagging
            (None, None) => Ok(Self::new(1, n, false, seed, n_jobs)),
            // Active (su)bagging
            (Some(k), _) => {
                let spec = subsample_size.unwrap_or(match replace {
                    false => Left(n),
                    true => Right(0.5),
                });
                let size = match spec {
                    Left(count) => count,
                    Right(share) => usize::max((share * n as f64).round() as usize, 1),
                };
                Ok(Self::new(k, size, replace, seed, n_jobs))
            }
        }
    }
}

pub struct Interpolation<'a, I: CovariateInterpolator + 'a> {
    interpolations: Vec<I>,
    threshold_map: &'a [ResponseCoordinate],
}

impl<I: CovariateInterpolator> CovariateInterpolator for Interpolation<'_, I> {
    fn interpolate_index(&self, index: usize) -> f64 {
        let n_subfits = self.interpolations.len();
        debug_assert_eq!(self.threshold_map.len() % n_subfits, 0);

        let subfit_thresholds = &self.threshold_map[index * n_subfits..];
        self.interpolations
            .iter()
            .zip(subfit_thresholds)
            .map(|(interpolation, &response)| interpolation.interpolate(response))
            .sum::<f64>()
            / n_subfits as f64
    }

    fn is_empty(&self) -> bool {
        self.threshold_map.is_empty()
    }

    fn len(&self) -> usize {
        let n_subfits = self.interpolations.len();
        debug_assert_eq!(self.threshold_map.len() % n_subfits, 0);
        // n_thresholds
        self.threshold_map.len() / n_subfits
    }
}

impl<I: CovariateInterpolator> IntoIterator for Interpolation<'_, I> {
    type Item = f64;
    type IntoIter = CdfInterpolation<Self>;

    fn into_iter(self) -> Self::IntoIter {
        CdfInterpolation::new(self)
    }
}

/// Granularity with which the outer progress tracker gets updated.
const TRACKER_STEPS: usize = 100;

struct MultiTracker<'a, const STEPS: usize> {
    parent: &'a dyn ProgressTracker,
    n_subsamples: usize,
    state: Mutex<State>,
}
struct State {
    /// For how many sub fits we know the exact size
    known: usize,
    /// Total number of tasks of all known subfits
    total: usize,
    /// Total number of tasks completed of all known subfits
    completed: usize,
    /// How many times we already incremented the parent
    ///
    /// This is needed because the prognosis of where we stand (the computed step) can jump forward
    /// or backward a bit when a new exact size becomes known.
    latest_step: usize,
}
impl State {
    /// What share of the total work has been done on a scale of `0..STEPS`?
    /// Like `(completed / total) * (known * n_subsamples) * STEPS`, then taking a floor.
    /// Returns 0 before any inner fit has reported its total.
    fn step<const STEPS: usize>(&self, n_subsamples: usize) -> usize {
        // can take large values, but should be fine with `usize` typically being 64 bits
        self.completed * self.known * STEPS / (self.total * n_subsamples)
    }
}
impl<'a, const STEPS: usize> MultiTracker<'a, STEPS> {
    fn new(n_subsamples: usize, parent: &'a dyn ProgressTracker) -> Self {
        parent.set_total(STEPS);
        Self {
            parent,
            n_subsamples,
            state: Mutex::new(State {
                known: 0,
                total: 0,
                completed: 0,
                latest_step: 0,
            }),
        }
    }
}
impl<const STEPS: usize> ProgressTracker for MultiTracker<'_, STEPS> {
    fn set_total(&self, total: usize) {
        let mut state_mut = self.state.lock().unwrap();
        state_mut.known += 1;
        state_mut.total += total;
    }

    fn increment(&self) {
        let (before, after) = {
            let mut state_mut = self.state.lock().unwrap();

            let before = state_mut.latest_step;
            state_mut.completed += 1;
            let after = state_mut.step::<STEPS>(self.n_subsamples);
            state_mut.latest_step = state_mut.latest_step.max(after);

            (before, after)
        };

        for _ in before..after {
            // this is typically 1 and never is a large amount
            self.parent.increment();
        }
    }

    /// irrelevant
    fn finish(&self) {}
}

macro_rules! impl_idr_fit_for {
    ($inner:ty) => {
        impl IsotonicDistributionalRegressionFit for Fit<$inner> {
            type Covariate<'a> = <$inner as IsotonicDistributionalRegressionFit>::Covariate<'a>;
            type CovariateOrder = <$inner as IsotonicDistributionalRegressionFit>::CovariateOrder;
            type Config = (
                Config,
                <$inner as IsotonicDistributionalRegressionFit>::Config,
            );

            fn fit(
                x: &[f64],
                y: &[f64],
                y_observed: Option<&[bool]>,
                weight: Option<&[f64]>,
                covariate_order: Self::CovariateOrder,
                response_order: StochasticOrder,
                decreasing: bool,
                config: Self::Config,
                progress: &dyn ProgressTracker,
            ) -> Result<Self, Error> {
                if y.len() == 0 || x.len() % y.len() != 0 {
                    return Err(Error::IncompatibleShapes {
                        covariate_len: x.len(),
                        response_len: y.len(),
                        weight_len: weight.map(|slice| slice.len()),
                        y_observed_len: y_observed.map(|slice| slice.len()),
                    });
                }
                let dimension = x.len() / y.len();
                let n = validate(x.chunks_exact(dimension), y, y_observed, weight)?;

                let (config, base_config) = config;

                let tracker = MultiTracker::<TRACKER_STEPS>::new(config.n_subsamples, progress);

                let fits = if config.subsample_size == n {
                    vec![<$inner>::fit(
                        x,
                        y,
                        y_observed,
                        weight,
                        covariate_order.clone(),
                        response_order,
                        decreasing,
                        base_config,
                        &tracker,
                    )?]
                } else {
                    // Generate all subsample index vectors serially for determinism.
                    fn gen_subsample<R: rand::Rng>(
                        rng: &mut R,
                        replace: bool,
                        n: usize,
                        subsample_size: usize,
                    ) -> Vec<usize> {
                        if replace {
                            let range = Uniform::new(0, n).unwrap();
                            (0..subsample_size).map(|_| range.sample(rng)).collect()
                        } else {
                            index::sample(rng, n, subsample_size.into()).into_vec()
                        }
                    }
                    let subsamples: Vec<Vec<usize>> = match config.seed {
                        Some(s) => {
                            use rand::SeedableRng;
                            let mut rng = rand::rngs::StdRng::seed_from_u64(s);
                            (0..config.n_subsamples)
                                .map(|_| {
                                    gen_subsample(
                                        &mut rng,
                                        config.replace,
                                        n,
                                        config.subsample_size,
                                    )
                                })
                                .collect()
                        }
                        None => {
                            let mut rng = rand::rng();
                            (0..config.n_subsamples)
                                .map(|_| {
                                    gen_subsample(
                                        &mut rng,
                                        config.replace,
                                        n,
                                        config.subsample_size,
                                    )
                                })
                                .collect()
                        }
                    };

                    let fit_subsample = |subsample: Vec<usize>| -> Result<$inner, Error> {
                        let len = config.subsample_size * dimension;
                        let mut sub_x = vec![f64::NAN; len];
                        let mut sub_y = vec![f64::NAN; config.subsample_size];
                        for (i, &sample) in subsample.iter().enumerate() {
                            sub_x[i * dimension..(i + 1) * dimension]
                                .copy_from_slice(&x[sample * dimension..(sample + 1) * dimension]);
                            sub_y[i] = y[sample];
                        }
                        assert!(sub_x.iter().all(|v| !v.is_nan()));
                        assert!(sub_y.iter().all(|v| !v.is_nan()));

                        let mut sub_y_observed = None;
                        if let Some(observed) = y_observed {
                            sub_y_observed =
                                Some(subsample.iter().map(|&i| observed[i]).collect::<Vec<_>>());
                        }

                        let mut sub_weight = None;
                        if let Some(weight) = weight {
                            sub_weight =
                                Some(subsample.iter().map(|&i| weight[i]).collect::<Vec<_>>());
                        }

                        <$inner>::fit(
                            &sub_x,
                            &sub_y,
                            sub_y_observed.as_deref(),
                            sub_weight.as_deref(),
                            covariate_order.clone(),
                            response_order,
                            decreasing,
                            base_config.clone(),
                            &tracker,
                        )
                    };

                    #[cfg(feature = "parallel")]
                    {
                        if config.n_jobs <= 1 {
                            subsamples
                                .into_iter()
                                .map(fit_subsample)
                                .collect::<Result<Vec<_>, _>>()?
                        } else {
                            let pool = rayon::ThreadPoolBuilder::new()
                                .num_threads(config.n_jobs)
                                .build()
                                .expect("failed to build rayon thread pool");
                            pool.install(|| {
                                subsamples
                                    .into_par_iter()
                                    .map(fit_subsample)
                                    .collect::<Result<Vec<_>, _>>()
                            })?
                        }
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        subsamples
                            .into_iter()
                            .map(fit_subsample)
                            .collect::<Result<Vec<_>, _>>()?
                    }
                };

                let covariates = unique_covariates(x, dimension);

                let mut thresholds = match y_observed {
                    Some(observed) => y
                        .iter()
                        .zip(observed)
                        .filter(|&(_, &o)| o)
                        .map(|(&r, _)| r)
                        .collect(),
                    None => Vec::from(y),
                };
                thresholds.sort_unstable_by(|a, b| a.total_cmp(b));
                thresholds.dedup();
                let threshold_map = derive_threshold_map(&thresholds, &fits);
                Ok(Self {
                    fits,
                    covariates,
                    thresholds,
                    threshold_map,
                    covariate_groups: covariate_order,
                })
            }

            fn interpolate_covariate<'a>(
                &'a self,
                covariate: Self::Covariate<'_>,
            ) -> impl CovariateInterpolator + IntoCdfIterator + 'a {
                let interpolations = self
                    .fits
                    .iter()
                    .map(|fit| fit.interpolate_covariate(covariate))
                    .collect();
                Interpolation {
                    interpolations,
                    threshold_map: &self.threshold_map,
                }
            }

            fn thresholds(&self) -> &[f64] {
                &self.thresholds
            }

            fn assert_consistent(&self) {
                for fit in &self.fits {
                    fit.assert_consistent();
                }
            }
        }
    };
}

#[cfg(feature = "partial-order")]
impl_idr_fit_for!(partial_order::Fit);
impl_idr_fit_for!(total_order::Fit);

#[cfg(feature = "partial-order")]
impl Fit<partial_order::Fit> {
    pub fn interpolate_covariate_with_workspace(
        &self,
        covariate: <Self as IsotonicDistributionalRegressionFit>::Covariate<'_>,
        workspace: &mut PredictionWorkspace,
    ) -> Interpolation<'_, partial_order::Interpolation<'_>> {
        let interpolations = self
            .fits
            .iter()
            .map(|fit| fit.interpolate_covariate_with_workspace(covariate, workspace))
            .collect();
        Interpolation {
            interpolations,
            threshold_map: &self.threshold_map,
        }
    }
}

impl Fit<total_order::Fit> {
    pub fn predict_grid<I: IntoIterator<Item = f64>>(
        &self,
        covariates: I,
        thresholds: impl IntoIterator<Item = f64>,
    ) -> GridPredictor<'_, I::IntoIter> {
        let mut covariate_iter = covariates.into_iter();
        let first_query = covariate_iter.next().expect("empty grid is not supported");
        let states = self
            .fits
            .iter()
            .map(|fit| {
                let mut search = CovariateSearch::new(&fit.covariates);
                let coordinate = search.advance(first_query);
                let cdfs = (fit.cdfs.as_slice(), fit.n_threshold());
                let interpolation = SingleInterpolation::from_coordinate(coordinate, cdfs);
                GridPredictorState {
                    search,
                    interpolation,
                    cdfs,
                }
            })
            .collect();
        let thresholds = search_responses_sorted(thresholds, &self.thresholds).collect();
        GridPredictor {
            covariate_iter,
            states,
            thresholds,
            threshold_index: 0,
            threshold_map: &self.threshold_map,
        }
    }
}

pub struct GridPredictor<'a, I> {
    covariate_iter: I,
    states: Vec<GridPredictorState<'a>>,
    thresholds: Vec<ResponseCoordinate>,
    threshold_index: usize,
    threshold_map: &'a [ResponseCoordinate],
}

impl<I: Iterator<Item = f64>> Iterator for GridPredictor<'_, I> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.threshold_index < self.thresholds.len() {
            let active_threshold = self.thresholds[self.threshold_index];
            let value = match active_threshold {
                ResponseCoordinate::StrictlyBelowAll => 0.0,
                ResponseCoordinate::AboveOrAtIndex(index) => {
                    let n_subagg = self.states.len();
                    self.states
                        .iter()
                        .zip(&self.threshold_map[index * n_subagg..])
                        .map(|(state, &response)| state.interpolation.interpolate(response))
                        .sum::<f64>()
                        / n_subagg as f64
                }
            };

            self.threshold_index += 1;
            if self.threshold_index == self.thresholds.len()
                && let Some(query) = self.covariate_iter.next()
            {
                for state in &mut self.states {
                    state.update(query);
                }
                self.threshold_index = 0;
            }

            Some(value)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (group_lower, group_upper) = self.covariate_iter.size_hint();
        let compute = |count| (count + 1) * self.thresholds.len() - self.threshold_index;
        let lower = compute(group_lower);
        let upper = group_upper.map(compute);
        (lower, upper)
    }
}

impl<I: ExactSizeIterator<Item = f64>> ExactSizeIterator for GridPredictor<'_, I> {}

fn unique_covariates(covariates: &[f64], dimension: usize) -> Vec<f64> {
    debug_assert_eq!(covariates.len() % dimension, 0);
    let n_covariates = covariates.len() / dimension;

    // Indices of rows: row i corresponds to covariates[i*n .. (i+1)*n]
    let mut rows = lexicographic_order(covariates, n_covariates, dimension);
    rows.dedup_by(|&mut i, &mut j| {
        let a = &covariates[i * dimension..(i + 1) * dimension];
        let b = &covariates[j * dimension..(j + 1) * dimension];
        a.iter().zip(b).all(|(aa, bb)| aa == bb)
    });

    let mut out = Vec::with_capacity(rows.len() * dimension);
    for i in rows {
        out.extend_from_slice(&covariates[i * dimension..(i + 1) * dimension]);
    }

    out
}

fn derive_threshold_map<F: IsotonicDistributionalRegressionFit>(
    thresholds: &[f64],
    fits: &[F],
) -> Vec<ResponseCoordinate> {
    // TODO: Avoid this transpose
    let mut indices: Vec<_> = fits
        .iter()
        .flat_map(|fit| search_responses_sorted(thresholds.iter().copied(), fit.thresholds()))
        .collect();
    transpose(&mut indices, fits.len(), thresholds.len());
    indices
}

#[cfg(test)]
mod test {
    use crate::partial_order::CovariateGroups;
    use crate::structures::{Increasing, StochasticOrder};
    use crate::subagging::{Config, Fit};
    use crate::{IsotonicDistributionalRegressionFit, partial_order, total_order};

    #[test]
    fn test_sd_symmetry_through_subagging() {
        let covariates = vec![
            0.2, 0.6, // row 0
            0.6, 0.2, // row 1  → same multiset as row 0
            0.9, 0.8, // row 2
            0.8, 0.9, // row 3  → same multiset as row 2
        ];
        let responses = vec![0.0, 1.0, 0.5, 0.6];
        let weight = vec![1.0; 4];
        let covariate_order = CovariateGroups::parse([("sd", [0usize, 1])], 2).unwrap();

        let fit = Fit::<partial_order::Fit>::fit(
            &covariates,
            &responses,
            None,
            Some(&weight),
            covariate_order,
            StochasticOrder::StochasticDominance,
            false,
            (Config::new(1, 4, false, None, 1), Default::default()),
            &crate::NoProgress,
        )
        .unwrap();

        // Check that SD-equivalent covariates get the same CDF
        let cdf_0_2_0_6: Vec<_> = fit.cdf(&[0.2, 0.6]).collect();
        let cdf_0_6_0_2: Vec<_> = fit.cdf(&[0.6, 0.2]).collect();
        assert_eq!(
            cdf_0_2_0_6, cdf_0_6_0_2,
            "SD-equivalent rows 0,1 must match"
        );

        let cdf_0_9_0_8: Vec<_> = fit.cdf(&[0.9, 0.8]).collect();
        let cdf_0_8_0_9: Vec<_> = fit.cdf(&[0.8, 0.9]).collect();
        assert_eq!(
            cdf_0_9_0_8, cdf_0_8_0_9,
            "SD-equivalent rows 2,3 must match"
        );
    }

    #[test]
    fn test_plain_wrapper() {
        let fit = Fit::<total_order::Fit>::fit(
            &[1.0, 2.0, 3.0],
            &[4.0, 3.0, 2.0],
            Some(&[false, true, true]),
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            (Config::new(1, 3, false, None, 1), Default::default()),
            &crate::NoProgress,
        )
        .unwrap();

        assert_eq!(fit.cdf(1.0).collect::<Vec<_>>(), vec![1.0 / 3.0, 2.0 / 3.0]);
        assert_eq!(fit.cdf(2.0).collect::<Vec<_>>(), vec![1.0 / 3.0, 2.0 / 3.0]);
        assert_eq!(fit.cdf(3.0).collect::<Vec<_>>(), vec![1.0 / 3.0, 2.0 / 3.0]);
        assert!(fit.mean(2.0).is_nan());
    }

    #[test]
    fn test_nans() {
        let covariates = [
            0.5488135, 0.71518937, 0.60276338, 0.54488318, 0.4236548, 0.64589411, 0.43758721,
            0.891773, 0.96366276, 0.38344152, 0.79172504, 0.52889492, 0.56804456, 0.92559664,
            0.07103606, 0.0871293, 0.0202184, 0.83261985, 0.77815675, 0.87001215, 0.97861834,
            0.79915856, 0.46147936, 0.78052918, 0.11827443, 0.63992102, 0.14335329, 0.94466892,
            0.52184832, 0.41466194, 0.26455561, 0.77423369, 0.45615033, 0.56843395, 0.0187898,
            0.6176355, 0.61209572, 0.616934, 0.94374808, 0.6818203, 0.3595079, 0.43703195,
            0.6976312, 0.06022547, 0.66676672, 0.67063787, 0.21038256, 0.1289263, 0.31542835,
            0.36371077, 0.57019677, 0.43860151, 0.98837384, 0.10204481, 0.20887676, 0.16130952,
            0.65310833, 0.2532916, 0.46631077, 0.24442559, 0.15896958, 0.11037514, 0.65632959,
            0.13818295, 0.19658236, 0.36872517, 0.82099323, 0.09710128, 0.83794491, 0.09609841,
            0.97645947, 0.4686512, 0.97676109, 0.60484552, 0.73926358, 0.03918779, 0.28280696,
            0.12019656, 0.2961402, 0.11872772, 0.31798318, 0.41426299, 0.0641475, 0.69247212,
            0.56660145, 0.26538949, 0.52324805, 0.09394051, 0.5759465, 0.9292962,
        ];
        let responses = [
            0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0,
            2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0,
        ];

        let fit = Fit::<partial_order::Fit>::fit(
            &covariates,
            &responses,
            None,
            None,
            CovariateGroups::empty(3),
            StochasticOrder::StochasticDominance,
            false,
            (Config::disable(responses.len()), Default::default()),
            &crate::NoProgress,
        )
        .unwrap();
        assert!(
            covariates
                .chunks_exact(3)
                .map(|c| fit.mean(c))
                .all(|v| !v.is_nan() && v.is_finite() && v >= 0.0 && v <= 2.0)
        );
    }

    #[test]
    fn test_larger() {
        let fit = Fit::<total_order::Fit>::fit(
            &(0..20).map(|i| i as f64).collect::<Vec<_>>(),
            &(10..30).map(|i| i as f64).collect::<Vec<_>>(),
            None,
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            (Config::new(5, 10, false, None, 1), Default::default()),
            &crate::NoProgress,
        )
        .unwrap();
        assert_eq!(fit.cdf_at(5.0, 5.0), 0.0);
        assert_eq!(fit.cdf_at(5.0, 35.0), 1.0);
        assert_eq!(fit.cdf_at(15.0, 5.0), 0.0);
        assert_eq!(fit.cdf_at(15.0, 35.0), 1.0);
    }

    #[test]
    fn test_small() {
        let covariates = [
            29.90845126,
            19.0484518,
            18.91235784,
            33.44648179,
            23.27694686,
        ];
        let responses = [3., 23., 165., 5., 57.];
        let fit = Fit::<total_order::Fit>::fit(
            &covariates,
            &responses,
            Some(&[false, false, true, false, true]),
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            true,
            (Config::new(1, 5, false, None, 1), Default::default()),
            &crate::NoProgress,
        )
        .unwrap();
        let result = covariates
            .into_iter()
            .flat_map(|c| fit.cdf(c))
            .collect::<Vec<_>>();
        assert_eq!(
            result,
            [
                1.0,
                1.0,
                (19.0484518 - 18.91235784) / (23.27694686 - 18.91235784),
                1.0,
                0.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ],
        );
    }
}
