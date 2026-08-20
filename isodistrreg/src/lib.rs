mod error;
mod float;
#[doc(hidden)]
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
pub use float::Float;
pub use prediction::{CovariateInterpolator, IntoCdfIterator, ResponseCoordinate, quantile};
pub use progress::{NoProgress, ProgressTracker};

/// Implementations of Isotonic Distributional Regression implement at least this functionality.
///
/// Whether based on covariates with a total order or partial order, these are the shared behaviors.
///
/// The associated types `X` and `Y` carry the precision of the user's covariate and response data.
/// Both are stored in the fit at that precision (no implicit f64 widening), and prediction inputs
/// likewise use the same precision. CDF *outputs* are always `f32` because the algorithm body
/// computes in `f32` post-preprocessing.
pub trait IsotonicDistributionalRegressionFit: Sized {
    /// Covariate scalar precision.
    type X: Float;
    /// Response/threshold scalar precision.
    type Y: Float;
    /// Prediction-input shape — `Self::X` for total-order fits, `&'a [Self::X]` for partial-order.
    type XInput<'a>;
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
    ///
    /// Weights `W` are independent of the covariate (`X`) and response (`Y`) precisions — the
    /// caller can pass `f32` or `f64`-typed weights regardless of `X`/`Y`. The implementation
    /// narrows once on read to whatever precision its algorithm operates in (`f32` for the
    /// total-order kernels, `f64` for the partial-order solver).
    ///
    /// Since the total-order kernels run their weight arithmetic in f32, it is up to the
    /// caller to keep the weights (and their sum) within f32's finite range and not too
    /// imbalanced — as a rule of thumb, no weight should be smaller than roughly 2⁻²⁴
    /// (~1e-7) of the total. Beyond that, small weights can lose their mass to round-off or
    /// the fit can produce NaN values. Positive weights that round to 0.0f32 are dropped
    /// like zero weights.
    ///
    /// Zero weights are allowed; such observations are dropped during preprocessing, so the
    /// fit equals the fit on the positive-weight subsample — including the threshold grid
    /// and covariate set. With every weight zero, the result is the empty fit (a sub-CDF
    /// that is 0 everywhere).
    #[allow(clippy::too_many_arguments)]
    fn fit<W: Float>(
        x: &[Self::X],
        y: &[Self::Y],
        y_observed: Option<&[bool]>,
        sample_weight: Option<&[W]>,
        x_order: Self::CovariateOrder,
        y_order: StochasticOrder,
        decreasing: bool,
        config: Self::Config,
        progress: &dyn ProgressTracker,
    ) -> Result<Self, Error>;

    fn interpolate_covariate<'a>(
        &'a self,
        x: Self::XInput<'_>,
    ) -> impl CovariateInterpolator + IntoCdfIterator + 'a;

    fn get_response_coordinate(&self, y: Self::Y) -> ResponseCoordinate {
        prediction::search_response(y, self.thresholds())
    }

    /// Predict the mean for a single covariate.
    ///
    /// If the fit is based on (partially) censored observations, the mean is not guaranteed to be
    /// finite. The integral is accumulated in `f64` internally for numerical stability (worst
    /// case n=100k thresholds, ~3 f32 ulps of running-sum round-off) and narrowed to `Self::Y`
    /// at the return.
    fn mean(&self, x: Self::XInput<'_>) -> Self::Y {
        prediction::mean(self.cdf(x), self.thresholds().iter().copied())
    }

    /// Predict the full (sub-)CDF at threshold in `thresholds()`.
    ///
    /// CDF values are returned as **f32** — that's the precision the algorithm body computes
    /// and stores. Widening here would pass through synthetic precision without buying
    /// anything.
    fn cdf(&self, x: Self::XInput<'_>) -> impl ExactSizeIterator<Item = f32> {
        self.interpolate_covariate(x).into_iter()
    }

    /// Predict the (sub-)CDF at specified threshold.
    fn cdf_at(&self, x: Self::XInput<'_>, y: Self::Y) -> f32 {
        let interpolation = self.interpolate_covariate(x);
        let y_coordinate = self.get_response_coordinate(y);
        interpolation.interpolate(y_coordinate)
    }

    fn quantile(&self, x: Self::XInput<'_>, probability: f32, upper: bool) -> Self::Y {
        let interpolator = self.interpolate_covariate(x);
        quantile(&interpolator, probability, upper, self.thresholds())
    }

    fn n_threshold(&self) -> usize {
        self.thresholds().len()
    }

    fn thresholds(&self) -> &[Self::Y];

    fn assert_consistent(&self);
}
