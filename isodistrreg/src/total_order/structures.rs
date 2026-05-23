use crate::Float;
use crate::error::Error;
use crate::prediction::{CovariateInterpolator, ResponseCoordinate, search_responses_sorted};
use crate::preprocessing::validate;
use crate::routines::transpose;
use crate::structures::{Increasing, Observation, StochasticOrder};
use crate::total_order::hazard_rate_order::preprocess_censored as preprocess_hazard_rate_censored;
use crate::total_order::prediction::{CovariateSearch, GridPredictorState, Interpolation};
use crate::total_order::preprocessing::preprocess_uncensored;
use crate::total_order::stochastic_dominance::preprocess_censored as preprocess_sd_censored;
use crate::total_order::weight_noise_floor;
use crate::total_order::{hazard_rate_order, stochastic_dominance};
use crate::{Decreasing, IntoCdfIterator, IsotonicDistributionalRegressionFit, ProgressTracker};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A computed IDR solution that can be used to make distributional predictions.
///
/// `X` and `Y` are the precision of covariate and response/threshold values respectively.
/// The CDF values are always stored in **f32**: the algorithm body (post-preprocessing) runs
/// in f32 throughout; widening to f64 here would waste memory and pass through synthetic
/// precision. See `total_order::stochastic_dominance::censored::fast` for the rationale.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Fit<X: Float, Y: Float> {
    /// Whether the fit is increasing (i.e., nondecreasing) w.r.t. the covariate ordering.
    pub increasing: bool,
    /// For each covariate, the estimated distribution.
    ///
    /// Represents a covariate-major matrix with in each minor index a CDF value for the
    /// corresponding threshold.
    pub cdfs: Vec<f32>,
    /// Unique covariate values sorted from low to high.
    pub covariates: Vec<X>,
    /// Unique thresholds sorted from low to high.
    pub thresholds: Vec<Y>,
    /// How well the numerical solve worked.
    pub quality_indicators: QualityIndicators,
}

/// Configuration for the total-order solver. Currently has no tunable knobs — the
/// f32 noise-clamping tolerance used by the hazard-rate and SD-censored kernels is
/// computed per-fit from the post-preprocessing observation count via
/// [`weight_noise_floor`].
#[derive(Clone, Default)]
pub struct Config;

impl<X: Float, Y: Float> IsotonicDistributionalRegressionFit for Fit<X, Y> {
    type X = X;
    type Y = Y;
    type XInput<'a> = X;
    /// The covariate order is increasing with increasing values.
    type CovariateOrder = Increasing;
    type Config = Config;

    fn fit<W: Float>(
        covariates: &[X],
        responses: &[Y],
        censoring: Option<&[bool]>,
        weights: Option<&[W]>,
        _covariate_order: Self::CovariateOrder,
        response_order: StochasticOrder,
        decreasing: bool,
        settings: Self::Config,
        progress: &dyn ProgressTracker,
    ) -> Result<Self, Error> {
        let n = validate(covariates.chunks_exact(1), responses, censoring, weights)?;

        let mut weights_allocation = None;
        let weights_to_use = weights.unwrap_or_else(|| {
            weights_allocation = Some(vec![W::one(); n]);
            weights_allocation.as_deref().unwrap()
        });

        // Preprocess once per case, then dispatch to the monomorphic f32 algorithm body. Each
        // algorithm runs against its pre-built `AlgorithmContext` (or `CensoredSdContext`); the
        // `W` narrowing happens entirely inside `preprocess_*`.
        let output: AlgorithmOutput<X, Y> = match censoring {
            Some(censoring) if censoring.iter().any(|&b| b) => match response_order {
                StochasticOrder::StochasticDominance => {
                    if decreasing {
                        // TODO: Running the censoring algorithm in reverse is not yet properly
                        //  supported, so we reverse manually here. The negation happens on the
                        //  raw `&[X]` before preprocessing.
                        let negated_covariates = covariates
                            .iter()
                            .map(|&v| X::zero() - v)
                            .collect::<Vec<_>>();
                        let context = preprocess_sd_censored(
                            &negated_covariates,
                            responses,
                            censoring,
                            weights_to_use,
                        )?;
                        let mut cdfs =
                            stochastic_dominance::censored::<Increasing, _, _>(&context, progress);
                        let mut unique_covariates = context.unique_covariates;
                        let thresholds = context.thresholds;
                        for c in unique_covariates.iter_mut() {
                            *c = X::zero() - *c;
                        }
                        unique_covariates.reverse();
                        debug_assert!(unique_covariates.windows(2).all(|w| {
                            w[0].partial_cmp(&w[1]).unwrap() != std::cmp::Ordering::Greater
                        }));
                        let n_covariate = unique_covariates.len();
                        for threshold in cdfs.chunks_exact_mut(n_covariate) {
                            threshold.reverse();
                        }
                        AlgorithmOutput {
                            cdfs,
                            unique_covariates,
                            thresholds,
                        }
                    } else {
                        let context = preprocess_sd_censored(
                            covariates,
                            responses,
                            censoring,
                            weights_to_use,
                        )?;
                        let cdfs =
                            stochastic_dominance::censored::<Increasing, _, _>(&context, progress);
                        AlgorithmOutput {
                            cdfs,
                            unique_covariates: context.unique_covariates,
                            thresholds: context.thresholds,
                        }
                    }
                }
                StochasticOrder::HazardRateOrder => {
                    let context = preprocess_hazard_rate_censored(
                        covariates,
                        responses,
                        censoring,
                        weights_to_use,
                    );
                    let cdfs = if decreasing {
                        hazard_rate_order::censored::<Decreasing, _, _>(&context, progress)
                    } else {
                        hazard_rate_order::censored::<Increasing, _, _>(&context, progress)
                    };
                    AlgorithmOutput {
                        cdfs,
                        unique_covariates: context.unique_covariates,
                        thresholds: context.unique_responses,
                    }
                }
            },
            _ => match response_order {
                StochasticOrder::StochasticDominance => {
                    let context = preprocess_uncensored(covariates, responses, weights_to_use);
                    let cdfs = if decreasing {
                        stochastic_dominance::uncensored::<Decreasing, _, _>(&context, progress)
                    } else {
                        stochastic_dominance::uncensored::<Increasing, _, _>(&context, progress)
                    };
                    AlgorithmOutput {
                        cdfs,
                        unique_covariates: context.unique_covariates,
                        thresholds: context.unique_responses,
                    }
                }
                StochasticOrder::HazardRateOrder => {
                    let context = preprocess_uncensored(covariates, responses, weights_to_use);
                    let cdfs = if decreasing {
                        hazard_rate_order::uncensored::<Decreasing, _, _>(&context, progress)
                    } else {
                        hazard_rate_order::uncensored::<Increasing, _, _>(&context, progress)
                    };
                    AlgorithmOutput {
                        cdfs,
                        unique_covariates: context.unique_covariates,
                        thresholds: context.unique_responses,
                    }
                }
            },
        };

        let AlgorithmOutput {
            mut cdfs,
            unique_covariates,
            thresholds,
        } = output;
        transpose(&mut cdfs, thresholds.len(), unique_covariates.len());

        let _ = settings;
        Ok(Self {
            increasing: !decreasing,
            cdfs,
            covariates: unique_covariates,
            thresholds,
            quality_indicators: QualityIndicators {
                epsilon: weight_noise_floor(n) as f64,
            },
        })
    }

    fn interpolate_covariate<'a>(
        &'a self,
        covariate: Self::XInput<'_>,
    ) -> impl CovariateInterpolator + IntoCdfIterator + 'a {
        if self.is_empty() {
            debug_assert!(
                self.cdfs.is_empty(),
                "we want to use this empty slice as a dummy value"
            );
            Interpolation::Exact { cdf: &self.cdfs }
        } else {
            Interpolation::new(
                covariate,
                &self.covariates,
                (&self.cdfs, self.n_threshold()),
            )
        }
    }

    fn thresholds(&self) -> &[Y] {
        &self.thresholds
    }

    fn assert_consistent(&self) {
        assert!(!self.covariates.is_empty());
        assert!(self.covariates.windows(2).all(|w| w[0] < w[1]));

        assert!(!self.thresholds.is_empty());
        assert!(self.thresholds.windows(2).all(|w| w[0] < w[1]));

        assert_eq!(
            self.cdfs.len(),
            self.covariates.len() * self.thresholds.len()
        );
        assert!(self.cdfs.iter().all(|v| (0.0..=1.0).contains(v)));
        assert!(
            self.cdfs
                .chunks_exact(self.thresholds.len())
                .all(|cdf| cdf.is_sorted())
        );
    }
}

pub struct GridPredictor<'a, X, I> {
    covariate_iter: I,
    state: GridPredictorState<'a, X>,
    thresholds: Vec<ResponseCoordinate>,
    threshold_index: usize,
}

impl<X: Float, I: Iterator<Item = X>> Iterator for GridPredictor<'_, X, I> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.threshold_index < self.thresholds.len() {
            let value = self
                .state
                .interpolation
                .interpolate(self.thresholds[self.threshold_index]);

            self.threshold_index += 1;
            if self.threshold_index == self.thresholds.len()
                && let Some(query) = self.covariate_iter.next()
            {
                self.state.update(query);
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

impl<X: Float, I: ExactSizeIterator<Item = X>> ExactSizeIterator for GridPredictor<'_, X, I> {}

impl<X: Float, Y: Float> Fit<X, Y> {
    /// An empty instance represents a sub-CDF that stays at value 0.0, always.
    fn is_empty(&self) -> bool {
        let result = self.cdfs.is_empty();
        debug_assert_eq!(self.covariates.is_empty(), result);
        debug_assert_eq!(self.thresholds.is_empty(), result);
        result
    }
}

impl<X: Float, Y: Float> Fit<X, Y> {
    pub fn predict_grid<I: IntoIterator<Item = X>>(
        &self,
        covariates: I,
        thresholds: impl IntoIterator<Item = Y>,
    ) -> GridPredictor<'_, X, I::IntoIter> {
        let mut search = CovariateSearch::new(&self.covariates);
        let mut covariate_iter = covariates.into_iter();
        let first = covariate_iter.next().expect("empty grid is not supported");
        let first_coordinate = search.advance(first);
        let cdfs = (self.cdfs.as_slice(), self.n_threshold());
        let thresholds = search_responses_sorted(thresholds, &self.thresholds).collect();
        GridPredictor {
            covariate_iter,
            state: GridPredictorState {
                search,
                interpolation: Interpolation::from_coordinate(first_coordinate, cdfs),
                cdfs,
            },
            thresholds,
            threshold_index: 0,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct AlgorithmOutput<X, Y> {
    /// Algorithm output CDF values, in f32. See `Fit::cdfs` for the precision rationale.
    pub cdfs: Vec<f32>,
    pub unique_covariates: Vec<X>,
    pub thresholds: Vec<Y>,
}

#[derive(Clone)]
pub struct Partition<W, V> {
    /// Right boundary of the partition, excluding. So partition is [previous, ..., last], index
    pub index: usize,
    pub weight: W,
    pub value: V,
}

/// Partition used in the post-preprocessing total-order algorithms (accelerated PAVA, hazard
/// rate, SD-censored fast path), which all run in f32. The generic `Partition<W, V>` is also
/// used directly by `tonic_regression` (parametric on `F: Float`) — see that module.
pub type WeightedPartition = Partition<f32, f32>;

/// Information that is aggregated per covariate, like the total weight of all observations for a
/// covariate.
///
/// Weights are stored as **f32** — the entire post-preprocessing algorithm runs in f32, and
/// this struct is part of that domain. Preprocessing (which produces these) keeps its
/// per-observation sums in the caller's input weight precision `W` and downcasts to f32 at
/// the boundary, so `W = f64` inputs get f64-quality accumulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CovariateStatistic {
    /// Total weight of this covariate, must be positive.
    pub weight: f32,
    /// Cumulative weight of this covariate and all that came before it (last covariate statistic
    /// then has the total weight of all observations).
    pub cumulative_weight: f32,
}

/// Input required to run the algorithm.
///
/// Observation weights are stored as **f32** to match the post-preprocessing algorithm's
/// precision; covariate values stay in their input precision `X` and response values in their
/// input precision `Y` (they're only used for sort/compare, no arithmetic).
#[derive(Clone, Debug, PartialEq)]
pub struct AlgorithmContext<X, Y, S> {
    /// Holds `n` items sorted by response from low to high, uncensored < censored, and covariate
    pub observations: Vec<Observation<usize, Y, S, f32>>,
    /// Holds information about each covariate
    pub covariate_statistics: Vec<CovariateStatistic>,
    /// Covariates
    pub unique_covariates: Vec<X>,
    /// Thresholds
    pub unique_responses: Vec<Y>,
}

pub fn allocate_and_sort<S, R, I: Into<Observation<f64, R, S>>>(
    data: impl Iterator<Item = I>,
) -> Vec<Observation<f64, R, S>> {
    let mut allocated: Vec<_> = data.map(I::into).collect();
    allocated.sort_by(|a, b| a.x.total_cmp(&b.x));
    allocated
}

#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct QualityIndicators {
    /// Noise floor used by the post-preprocessing f32 kernels at fit time. Set by
    /// `weight_noise_floor(n_total)` during `Fit::fit`; for fits reconstructed via
    /// `from_cdfs` (where the original `n_total` is unknown) this is `f64::NAN`.
    pub epsilon: f64,
}

#[cfg(test)]
mod test {
    use crate::IsotonicDistributionalRegressionFit;
    use crate::structures::{Increasing, StochasticOrder};
    use crate::total_order::{Config, Fit, QualityIndicators};

    #[test]
    fn test_predict() {
        let fit: Fit<f64, f64> = Fit {
            increasing: true,
            cdfs: vec![
                0.25, 0.25, 0.5, 1.0, 0.25, 0.25, 0.25, 1.0, 0.0, 0.0, 0.25, 0.5,
            ],
            covariates: vec![0.0, 1.0, 2.0],
            thresholds: vec![0.0, 1.0, 2.0, 3.0],
            quality_indicators: QualityIndicators { epsilon: 0.0 },
        };

        // exact
        assert_eq!(
            fit.predict_grid([0., 1.], fit.thresholds.iter().copied())
                .collect::<Vec<_>>(),
            vec![0.25, 0.25, 0.5, 1.0, 0.25, 0.25, 0.25, 1.0],
        );
        assert_eq!(
            fit.predict_grid([0.0, 1.0], [-1.0, 1.0, 2.5, 3.0, 3.5])
                .collect::<Vec<_>>(),
            vec![0.0, 0.25, 0.5, 1.0, 1.0, 0.0, 0.25, 0.25, 1.0, 1.0],
        );

        // interpolating
        assert_eq!(
            fit.cdf(0.5).collect::<Vec<_>>(),
            vec![0.25, 0.25, 0.75 / 2.0, 1.0],
        );
        assert_eq!(fit.cdf_at(0.5, -1.0), 0.0);
        assert_eq!(fit.cdf_at(0.5, 1.0), 0.25);
        assert_eq!(fit.cdf_at(0.5, 2.5), 0.75 / 2.0);
        assert_eq!(fit.cdf_at(0.5, 3.0), 1.0);
        assert_eq!(fit.cdf_at(0.5, 3.5), 1.0);

        // out of range
        assert_eq!(
            fit.predict_grid([-0.5, 2.5], [0.0, 1.0, 2.0, 3.0])
                .collect::<Vec<_>>(),
            vec![0.25, 0.25, 0.5, 1.0, 0.0, 0.0, 0.25, 0.5],
        );
        assert_eq!(
            fit.predict_grid([-0.5, 2.5], [-1.0, 1.0, 2.5, 3.0, 3.5])
                .collect::<Vec<_>>(),
            vec![0.0, 0.25, 0.5, 1.0, 1.0, 0.0, 0.0, 0.25, 0.5, 0.5],
        );
    }

    /// All four (X, Y) precision combinations fit the same data to ~the same CDFs.
    ///
    /// `Y=f32` loses precision in thresholds (n ≤ 8 here, all integers, so no loss in practice);
    /// `X=f32` loses precision in the covariate grid (same). With this trivial test data the
    /// answers should be bit-equal for f32-vs-f64 covariates and tied to ~ulp for f32-vs-f64
    /// thresholds.
    #[test]
    fn test_all_four_xy_combinations() {
        let x_f64 = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y_f64 = [10.0, 30.0, 20.0, 50.0, 40.0];
        let x_f32: Vec<f32> = x_f64.iter().map(|&v| v as f32).collect();
        let y_f32: Vec<f32> = y_f64.iter().map(|&v| v as f32).collect();

        let fit_f64_f64 = Fit::<f64, f64>::fit::<f64>(
            &x_f64,
            &y_f64,
            None,
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();
        let fit_f32_f64 = Fit::<f32, f64>::fit::<f64>(
            &x_f32,
            &y_f64,
            None,
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();
        let fit_f64_f32 = Fit::<f64, f32>::fit::<f64>(
            &x_f64,
            &y_f32,
            None,
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();
        let fit_f32_f32 = Fit::<f32, f32>::fit::<f64>(
            &x_f32,
            &y_f32,
            None,
            None,
            Increasing,
            StochasticOrder::StochasticDominance,
            false,
            Config::default(),
            &crate::NoProgress,
        )
        .unwrap();

        // Cdfs should match across all four — integer covariates/thresholds round-trip exactly.
        assert_eq!(fit_f64_f64.cdfs, fit_f32_f64.cdfs);
        assert_eq!(fit_f64_f64.cdfs, fit_f64_f32.cdfs);
        assert_eq!(fit_f64_f64.cdfs, fit_f32_f32.cdfs);

        // Storage precision is preserved.
        assert_eq!(fit_f64_f64.covariates.len(), 5);
        assert_eq!(fit_f32_f64.covariates.len(), 5);
        assert_eq!(fit_f64_f32.thresholds.len(), 5);
        assert_eq!(fit_f32_f32.thresholds.len(), 5);
    }
}
