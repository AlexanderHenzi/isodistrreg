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
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "T: Serialize, T::X: Serialize, T::Y: Serialize, T::CovariateOrder: Serialize",
        deserialize = "T: Deserialize<'de>, T::X: Deserialize<'de>, T::Y: Deserialize<'de>, T::CovariateOrder: Deserialize<'de>"
    ))
)]
pub struct Fit<T: IsotonicDistributionalRegressionFit> {
    pub fits: Vec<T>,
    /// All covariates that were available to this problem, deduplicated.
    pub covariates: Vec<T::X>,
    /// All uncensored responses that were available to this problem, deduplicated.
    pub thresholds: Vec<T::Y>,
    /// For each fit, map the global threshold to the fit of the threshold.
    pub threshold_map: Vec<ResponseCoordinate>,
    pub covariate_groups: T::CovariateOrder,
}

impl<T: IsotonicDistributionalRegressionFit> Fit<T> {
    /// Creates a [`Fit`] directly from pre-computed components.
    pub fn from_parts(
        fit: T,
        covariates: Vec<T::X>,
        thresholds: Vec<T::Y>,
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
        debug_assert_ne!(n_subsamples, 0);
        debug_assert_ne!(subsample_size, 0);
        debug_assert_ne!(n_jobs, 0);

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

    /// Re-checks against the training-data size `n` the invariants that
    /// [`Config::parse`] enforces but the infallible [`Config::new`] cannot.
    ///
    /// Without this, a zero subsample count yields a fit whose every prediction
    /// divides by the member count, a subsample size above `n` panics inside
    /// `rand::seq::index::sample` when drawing without replacement, and a
    /// subsample size of zero surfaces as a confusing shape error from the
    /// inner fits.
    fn validate(&self, n: usize) -> Result<(), Error> {
        use crate::Error::SubaggingParameterInconsistency as BadParam;

        if self.n_subsamples == 0 {
            return Err(BadParam("subsample count should be >= 1"));
        }
        if self.subsample_size == 0 || self.subsample_size > n {
            return Err(BadParam(
                "subsample size should be at least 1 and at most the number of observations",
            ));
        }
        if self.n_jobs == 0 {
            return Err(BadParam("n_jobs must be at least 1"));
        }
        Ok(())
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

        fn count_from_share(share: f64, n: usize) -> usize {
            usize::max((share * n as f64).ceil() as usize, 1)
        }

        // A full-size draw WITHOUT replacement is a permutation of the training data, so
        // every such "subsample" reproduces the plain fit. (A full-size draw WITH
        // replacement is the classic bootstrap and is meaningful.)
        let full_without_replacement = !replace
            && match subsample_size {
                Some(Left(count)) => count == n,
                Some(Right(share)) => count_from_share(share, n) == n,
                None => false,
            };

        match (subsamples, subsample_size) {
            // Explicitly requesting the full sample without replacement is equivalent to
            // no subsampling at all; aggregating several identical fits would only waste
            // compute, so more than one subsample is rejected rather than silently
            // collapsed.
            (None | Some(1), Some(_)) if full_without_replacement => {
                Ok(Self::new(1, n, false, seed, n_jobs))
            }
            (Some(_), Some(_)) if full_without_replacement => Err(BadParam(
                "a full-size subsample without replacement reproduces the plain fit; \
                 choose a smaller subsample_size or set replace to sample with replacement",
            )),
            (None, Some(_)) => Err(BadParam(
                "specify the number of subsamples when specifying a subsample size",
            )),
            (None, None) if replace => Err(BadParam(
                "can't sample with replacement if subsample_size not specified",
            )),
            // No (su)bagging
            (None, None) => Ok(Self::new(1, n, false, seed, n_jobs)),
            // Active (su)bagging. Defaults: half the data for subagging (without
            // replacement), the full sample size for bootstrapping (with replacement).
            (Some(k), spec) => {
                let spec = spec.unwrap_or(match replace {
                    false => Right(0.5),
                    true => Right(1.0),
                });
                let size = match spec {
                    Left(count) => count,
                    Right(share) => count_from_share(share, n),
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
    fn interpolate_index(&self, index: usize) -> f32 {
        let n_subfits = self.interpolations.len();
        debug_assert_eq!(self.threshold_map.len() % n_subfits, 0);

        let subfit_thresholds = &self.threshold_map[index * n_subfits..];
        let sum: f32 = self
            .interpolations
            .iter()
            .zip(subfit_thresholds)
            .map(|(interpolation, &response)| interpolation.interpolate(response))
            .sum();
        sum / n_subfits as f32
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
    type Item = f32;
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
        impl<X: crate::Float, Y: crate::Float> IsotonicDistributionalRegressionFit for Fit<$inner>
        where
            <$inner as IsotonicDistributionalRegressionFit>::Config: Clone,
            <$inner as IsotonicDistributionalRegressionFit>::CovariateOrder: Clone,
        {
            type X = X;
            type Y = Y;
            type XInput<'a> = <$inner as IsotonicDistributionalRegressionFit>::XInput<'a>;
            type CovariateOrder = <$inner as IsotonicDistributionalRegressionFit>::CovariateOrder;
            type Config = (
                Config,
                <$inner as IsotonicDistributionalRegressionFit>::Config,
            );

            fn fit<W: crate::Float>(
                x: &[X],
                y: &[Y],
                y_observed: Option<&[bool]>,
                weight: Option<&[W]>,
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

                // `Config::new` is infallible and bypasses `Config::parse`'s validation;
                // an inconsistent configuration would otherwise panic deep in the
                // resampling code or yield a fit that panics on prediction.
                config.validate(n)?;

                let tracker = MultiTracker::<TRACKER_STEPS>::new(config.n_subsamples, progress);

                // A full-size subsample is the whole data set only when drawn WITHOUT
                // replacement; a with-replacement draw of size n is a genuine bootstrap
                // resample (it almost surely contains duplicates) and must go through
                // the resampling path below.
                let fits = if config.subsample_size == n && !config.replace {
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
                        let mut sub_x: Vec<X> = vec![X::nan(); len];
                        let mut sub_y: Vec<Y> = vec![Y::nan(); config.subsample_size];
                        for (i, &sample) in subsample.iter().enumerate() {
                            sub_x[i * dimension..(i + 1) * dimension]
                                .copy_from_slice(&x[sample * dimension..(sample + 1) * dimension]);
                            sub_y[i] = y[sample];
                        }
                        // Internal invariant (the copy loop filled every slot of the
                        // NaN-seeded buffers) — debug-only: the scan is O(subsample ·
                        // dimension) per subsample, and user-facing NaN rejection already
                        // happened in `validate` at the fit boundary.
                        debug_assert!(sub_x.iter().all(|v| !v.is_nan()));
                        debug_assert!(sub_y.iter().all(|v| !v.is_nan()));

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

                // Per the fit() contract zero-weight observations are inert: the
                // aggregate's covariate set and threshold grid must equal those of the
                // positive-weight subsample (the inner fits already guarantee this for
                // their own grids).
                let covariates = match weight {
                    None => unique_covariates(x, dimension),
                    Some(w) => {
                        let filtered: Vec<X> = (0..n)
                            .filter(|&i| w[i] > W::zero())
                            .flat_map(|i| x[i * dimension..(i + 1) * dimension].iter().copied())
                            .collect();
                        unique_covariates(&filtered, dimension)
                    }
                };

                let mut thresholds: Vec<Y> = (0..n)
                    .filter(|&i| weight.is_none_or(|w| w[i] > W::zero()))
                    .filter(|&i| y_observed.is_none_or(|o| o[i]))
                    .map(|i| y[i])
                    .collect();
                thresholds.sort_unstable_by(|a, b| crate::functionals::TotalCmp::total_cmp(a, b));
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
                covariate: Self::XInput<'_>,
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

            fn thresholds(&self) -> &[Y] {
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
impl_idr_fit_for!(partial_order::Fit<X, Y>);
impl_idr_fit_for!(total_order::Fit<X, Y>);

#[cfg(feature = "partial-order")]
impl<X: crate::Float, Y: crate::Float> Fit<partial_order::Fit<X, Y>> {
    pub fn interpolate_covariate_with_workspace(
        &self,
        covariate: <Self as IsotonicDistributionalRegressionFit>::XInput<'_>,
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

impl<X: crate::Float, Y: crate::Float> Fit<total_order::Fit<X, Y>> {
    pub fn predict_grid<I: IntoIterator<Item = X>>(
        &self,
        covariates: I,
        thresholds: impl IntoIterator<Item = Y>,
    ) -> GridPredictor<'_, X, I::IntoIter> {
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

pub struct GridPredictor<'a, X, I> {
    covariate_iter: I,
    states: Vec<GridPredictorState<'a, X>>,
    thresholds: Vec<ResponseCoordinate>,
    threshold_index: usize,
    threshold_map: &'a [ResponseCoordinate],
}

impl<X: crate::Float, I: Iterator<Item = X>> Iterator for GridPredictor<'_, X, I> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.threshold_index < self.thresholds.len() {
            let active_threshold = self.thresholds[self.threshold_index];
            let value: f32 = match active_threshold {
                ResponseCoordinate::StrictlyBelowAll => 0.0,
                ResponseCoordinate::AboveOrAtIndex(index) => {
                    let n_subagg = self.states.len();
                    let sum: f32 = self
                        .states
                        .iter()
                        .zip(&self.threshold_map[index * n_subagg..])
                        .map(|(state, &response)| state.interpolation.interpolate(response))
                        .sum();
                    sum / n_subagg as f32
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

impl<X: crate::Float, I: ExactSizeIterator<Item = X>> ExactSizeIterator
    for GridPredictor<'_, X, I>
{
}

fn unique_covariates<X: crate::Float>(covariates: &[X], dimension: usize) -> Vec<X> {
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
    thresholds: &[F::Y],
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
    use itertools::Either;

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

        let fit = Fit::<partial_order::Fit<f64, f64>>::fit::<f64>(
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
        let fit = Fit::<total_order::Fit<f64, f64>>::fit::<f64>(
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

        let fit = Fit::<partial_order::Fit<f64, f64>>::fit::<f64>(
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
        let fit = Fit::<total_order::Fit<f64, f64>>::fit::<f64>(
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

    /// A subsample fraction of 1.0 combined with `replace = true` is the
    /// classic bootstrap and must keep the full subsample size `n`.
    #[test]
    fn parse_full_size_bootstrap_keeps_size() {
        let config = Config::parse(Some(7), Some(Either::Right(1.0)), true, 10, Some(3), 1)
            .expect("b=7, p=1.0, replace=true is a valid bootstrap configuration");
        assert_eq!(config.n_subsamples, 7);
        assert!(config.replace);
        assert_eq!(
            config.subsample_size, 10,
            "p = 1.0 must give subsamples of size n = 10 (got {})",
            config.subsample_size
        );
    }

    /// Fractional sizes resolve to ceiling(n * p), per the documented `idrbag`
    /// contract ("a random subsample of size ceiling(nrow(X)*p)").
    #[test]
    fn parse_fractional_size_uses_ceiling() {
        let config = Config::parse(Some(2), Some(Either::Right(0.61)), false, 5, None, 1).unwrap();
        assert_eq!(config.subsample_size, 4); // ceiling(3.05)
    }

    /// The aggregate's threshold grid and covariate set come from the positive-weight
    /// subsample, matching the member fits' own grids; with every weight zero the
    /// aggregate is the empty fit.
    #[test]
    fn zero_weight_observations_are_inert_in_aggregate_grids() {
        let fit = |x: &[f64], y: &[f64], w: Option<&[f64]>| {
            Fit::<total_order::Fit<f64, f64>>::fit::<f64>(
                x,
                y,
                None,
                w,
                Increasing,
                StochasticOrder::StochasticDominance,
                false,
                (Config::disable(x.len()), Default::default()),
                &crate::NoProgress,
            )
            .unwrap()
        };
        let base = fit(&[1.0, 2.0], &[1.0, 2.0], None);
        let padded = fit(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0], Some(&[1.0, 1.0, 0.0]));
        assert_eq!(padded.thresholds(), base.thresholds());
        assert_eq!(padded.covariates, base.covariates);

        let empty = fit(&[1.0, 2.0], &[1.0, 2.0], Some(&[0.0, 0.0]));
        assert!(empty.thresholds().is_empty());
        assert_eq!(empty.cdf(1.0).count(), 0);
    }

    /// Without an explicit size, subagging (no replacement) defaults to half the data
    /// and bootstrapping (with replacement) to the full sample size.
    #[test]
    fn parse_default_sizes() {
        let subagged = Config::parse(Some(7), None, false, 10, None, 1).unwrap();
        assert_eq!(subagged.n_subsamples, 7);
        assert!(!subagged.replace);
        assert_eq!(subagged.subsample_size, 5);

        let bootstrap = Config::parse(Some(7), None, true, 10, None, 1).unwrap();
        assert_eq!(bootstrap.n_subsamples, 7);
        assert!(bootstrap.replace);
        assert_eq!(bootstrap.subsample_size, 10);
    }

    /// A full-size subsample without replacement is a permutation of the data: every
    /// subfit equals the plain fit, so aggregating more than one is rejected, while a
    /// single one collapses to the plain configuration.
    #[test]
    fn parse_full_size_without_replacement() {
        for spec in [Either::Left(10), Either::Right(1.0)] {
            assert!(Config::parse(Some(7), Some(spec), false, 10, None, 1).is_err());

            let single = Config::parse(Some(1), Some(spec), false, 10, None, 1).unwrap();
            assert_eq!(single.n_subsamples, 1);
            assert_eq!(single.subsample_size, 10);
            assert!(!single.replace);

            // Without a subsample count this remains the documented "no subagging" form.
            let plain = Config::parse(None, Some(spec), false, 10, None, 1).unwrap();
            assert_eq!(plain.n_subsamples, 1);
            assert_eq!(plain.subsample_size, 10);
        }
    }

    /// `Config::new` bypasses `parse`'s validation, so `fit` re-checks every
    /// invariant against the training-data size.
    #[test]
    fn validate_checks_all_invariants() {
        let config = |n_subsamples, subsample_size, replace, n_jobs| Config {
            n_subsamples,
            subsample_size,
            replace,
            seed: None,
            n_jobs,
        };

        assert!(config(2, 3, false, 1).validate(5).is_ok());
        // Classic bootstrap: full-size draws with replacement are valid.
        assert!(config(2, 5, true, 1).validate(5).is_ok());

        // Zero subsamples would yield a fit that divides by the member count.
        assert!(config(0, 3, false, 1).validate(5).is_err());
        // A subsample needs at least one observation.
        assert!(config(2, 0, false, 1).validate(5).is_err());
        // More draws than observations panic in the without-replacement sampler...
        assert!(config(2, 6, false, 1).validate(5).is_err());
        // ...and exceed the size contract of `parse` even with replacement.
        assert!(config(2, 6, true, 1).validate(5).is_err());
        // At least one worker thread is required.
        assert!(config(2, 3, false, 0).validate(5).is_err());
    }

    /// The fit-time consistency check turns the oversized-subsample panic inside
    /// `rand::seq::index::sample` into a recoverable error.
    #[test]
    fn fit_rejects_oversized_subsample() {
        let result = Fit::<total_order::Fit<f64, f64>>::fit::<f64>(
            &[1.0, 2.0, 3.0],
            &[4.0, 3.0, 2.0],
            None,
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            (Config::new(2, 4, false, Some(0), 1), Default::default()),
            &crate::NoProgress,
        );
        assert!(matches!(
            result,
            Err(crate::Error::SubaggingParameterInconsistency(_))
        ));
    }

    /// An all-censored subsample yields an empty member fit; the grid prediction
    /// path handles it like the pointwise path does (an empty member contributes 0
    /// everywhere).
    #[test]
    fn predict_grid_handles_empty_member_fits() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let y = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        // One event, five censored observations: a size-2 subsample without the
        // event produces the empty member fit.
        let observed = [true, false, false, false, false, false];
        let fit_with_seed = |seed| {
            Fit::<total_order::Fit<f64, f64>>::fit::<f64>(
                &x,
                &y,
                Some(&observed),
                None,
                Increasing,
                StochasticOrder::StochasticDominance,
                false,
                (Config::new(3, 2, false, Some(seed), 1), Default::default()),
                &crate::NoProgress,
            )
            .unwrap()
        };
        let agg = (0..1000)
            .map(fit_with_seed)
            .find(|agg| {
                agg.fits.iter().any(|f| f.thresholds.is_empty())
                    && agg.fits.iter().any(|f| !f.thresholds.is_empty())
            })
            .expect("some seed in 0..1000 must draw both member kinds");

        let via_cdf: Vec<f32> = x.iter().flat_map(|&c| agg.cdf(c)).collect();
        let via_grid: Vec<f32> = agg
            .predict_grid(x.iter().copied(), agg.thresholds().iter().copied())
            .collect();
        assert_eq!(via_grid, via_cdf);
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
        let fit = Fit::<total_order::Fit<f64, f64>>::fit::<f64>(
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
        let result: Vec<f32> = covariates.into_iter().flat_map(|c| fit.cdf(c)).collect();
        // Match the production path: prediction.rs computes the interpolation share in `X`
        // precision (here f64) and only narrows to f32 at the end.
        let share =
            ((19.0484518_f64 - 18.91235784_f64) / (23.27694686_f64 - 18.91235784_f64)) as f32;
        let expected: Vec<f32> = vec![1.0, 1.0, share, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        assert!(crate::test::is_relative_eq_vec(&result, &expected));
    }
}
