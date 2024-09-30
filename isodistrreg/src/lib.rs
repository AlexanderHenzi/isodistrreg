mod error;
pub mod functionals;
#[cfg(feature = "partial-order")]
pub mod partial_order;
mod prediction;
mod preprocessing;
pub mod progress;
pub mod routines;
mod structures;
#[cfg(feature = "subagging")]
pub mod subagging;
#[cfg(test)]
mod test;
pub mod total_order;

pub use crate::structures::{
    Decreasing, Direction, Increasing, Observation, Parallel, Serial, StochasticOrder,
};
pub use error::Error;
pub use prediction::{CovariateInterpolator, IntoCdfIterator, ResponseCoordinate, quantile};
pub use progress::{NoProgress, ProgressTracker};

/// Implementations of Isotonic Distributional Regression implement at least this functionality.
///
/// Whether based on covariates with a total order or partial order, these are the shared behaviors.
pub trait IsotonicDistributionalRegressionFit: Sized {
    type Covariate<'a>;
    type CovariateOrder: ?Sized;
    type Config;

    /// Fit an isotonic distributional regression.
    ///
    /// May take a long time to run depending on the input size, whether any inputs are censored and
    /// how many partitions the solution will need (i.e., how steadily in- / decreasing the response
    /// is with respect to the covariate).
    ///
    /// `progress` is called once per finished threshold during the fit. Pass [`NoProgress`] to
    /// disable progress reporting (the calls become indirect no-ops; overhead is negligible
    /// compared to per-threshold work).
    #[allow(clippy::too_many_arguments)]
    fn fit(
        x: &[f64],
        y: &[f64],
        y_observed: Option<&[bool]>,
        sample_weight: Option<&[f64]>,
        x_order: Self::CovariateOrder,
        y_order: StochasticOrder,
        decreasing: bool,
        config: Self::Config,
        progress: &dyn ProgressTracker,
    ) -> Result<Self, Error>;

    fn interpolate_covariate<'a>(
        &'a self,
        x: Self::Covariate<'_>,
    ) -> impl CovariateInterpolator + IntoCdfIterator + 'a;

    fn get_response_coordinate(&self, y: f64) -> ResponseCoordinate {
        prediction::search_response(y, self.thresholds())
    }

    /// Predict the mean for a single covariate.
    ///
    /// If the fit is based on (partially) censored observations, the mean is not guaranteed to be
    /// finite.
    fn mean(&self, x: Self::Covariate<'_>) -> f64 {
        prediction::mean(self.cdf(x), self.thresholds().iter().copied())
    }

    /// Predict the full (sub-)CDF at threshold in `thresholds()`.
    fn cdf(&self, x: Self::Covariate<'_>) -> impl ExactSizeIterator<Item = f64> {
        self.interpolate_covariate(x).into_iter()
    }

    /// Predict the (sub-)CDF at specified threshold.
    fn cdf_at(&self, x: Self::Covariate<'_>, y: f64) -> f64 {
        let interpolation = self.interpolate_covariate(x);
        let y_coordinate = self.get_response_coordinate(y);
        interpolation.interpolate(y_coordinate)
    }

    fn quantile(&self, x: Self::Covariate<'_>, probability: f64, upper: bool) -> f64 {
        let interpolator = self.interpolate_covariate(x);
        quantile(&interpolator, probability, upper, self.thresholds())
    }

    fn n_threshold(&self) -> usize {
        self.thresholds().len()
    }

    fn thresholds(&self) -> &[f64];

    fn assert_consistent(&self);
}
