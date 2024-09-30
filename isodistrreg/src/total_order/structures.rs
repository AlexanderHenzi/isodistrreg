use crate::error::Error;
use crate::prediction::{CovariateInterpolator, ResponseCoordinate, search_responses_sorted};
use crate::preprocessing::validate;
use crate::routines::transpose;
use crate::structures::{Increasing, Observation, StochasticOrder};
use crate::total_order::prediction::{CovariateSearch, GridPredictorState, Interpolation};
use crate::total_order::{hazard_rate_order, stochastic_dominance};
use crate::{Decreasing, IntoCdfIterator, IsotonicDistributionalRegressionFit, ProgressTracker};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A computed IDR solution that can be used to make distributional predictions.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Fit {
    /// Whether the fit is increasing (i.e., nondecreasing) w.r.t. the covariate ordering.
    pub increasing: bool,
    /// For each covariate, the estimated distribution.
    ///
    /// Represents a covariate-major matrix with in each minor index a CDF value for the
    /// corresponding threshold.
    pub cdfs: Vec<f64>,
    /// Unique covariate values sorted from low to high.
    pub covariates: Vec<f64>,
    /// Unique thresholds sorted from low to high.
    pub thresholds: Vec<f64>,
    /// How well the numerical solve worked.
    pub quality_indicators: QualityIndicators,
}

/// See the `Config` struct in the `any_order` module.
#[derive(Clone)]
pub struct Config {
    /// CDF values below `epsilon` or above `1 - epsilon` might be set to `0.0` or `1.0`
    /// respectively, used only to deal with numerical errors arising from operations on observation
    /// weights (like computing which share of weight has already been processed). Passing more
    /// balanced and normalized weights as input can allow this value to be lowered.
    pub epsilon: f64,
}
impl Default for Config {
    fn default() -> Self {
        Self { epsilon: 1e-8 }
    }
}

impl IsotonicDistributionalRegressionFit for Fit {
    type Covariate<'a> = f64;
    /// The covariate order is increasing with increasing values.
    type CovariateOrder = Increasing;
    type Config = Config;

    fn fit(
        covariates: &[f64],
        responses: &[f64],
        censoring: Option<&[bool]>,
        weights: Option<&[f64]>,
        _covariate_order: Self::CovariateOrder,
        response_order: StochasticOrder,
        decreasing: bool,
        settings: Self::Config,
        progress: &dyn ProgressTracker,
    ) -> Result<Self, Error> {
        let n = validate(covariates.chunks_exact(1), responses, censoring, weights)?;

        let mut weights_allocation = None;
        let weights_to_use = weights.unwrap_or_else(|| {
            weights_allocation = Some(vec![1.0; n]);
            weights_allocation.as_deref().unwrap()
        });

        // The different cases below use different pre-processing, so it happens in the different
        // `algorithm` functions.

        let output = match censoring {
            Some(censoring) if censoring.iter().any(|&b| b) => match response_order {
                StochasticOrder::StochasticDominance => {
                    if decreasing {
                        // TODO: Running the censoring algorithm in reverse is not yet properly
                        //  supported, so we reverse manually here
                        let negated_covariates = covariates.iter().map(|&v| -v).collect::<Vec<_>>();
                        let mut output = stochastic_dominance::censored::<Increasing>(
                            &negated_covariates,
                            responses,
                            censoring,
                            weights_to_use,
                            progress,
                        );

                        if let Ok(AlgorithmOutput {
                            cdfs,
                            unique_covariates,
                            ..
                        }) = &mut output
                        {
                            for c in unique_covariates.iter_mut() {
                                *c = -*c;
                            }
                            unique_covariates.reverse();
                            debug_assert!(unique_covariates.is_sorted());
                            let n_covariate = unique_covariates.len();
                            for threshold in cdfs.chunks_exact_mut(n_covariate) {
                                threshold.reverse();
                            }
                        }

                        output
                    } else {
                        stochastic_dominance::censored::<Increasing>(
                            covariates,
                            responses,
                            censoring,
                            weights_to_use,
                            progress,
                        )
                    }
                }
                StochasticOrder::HazardRateOrder => {
                    if decreasing {
                        hazard_rate_order::censored::<Decreasing>(
                            covariates,
                            responses,
                            censoring,
                            weights_to_use,
                            settings.epsilon,
                            progress,
                        )
                    } else {
                        hazard_rate_order::censored::<Increasing>(
                            covariates,
                            responses,
                            censoring,
                            weights_to_use,
                            settings.epsilon,
                            progress,
                        )
                    }
                }
            },
            _ => match response_order {
                StochasticOrder::StochasticDominance => {
                    if decreasing {
                        stochastic_dominance::uncensored::<Decreasing>(
                            covariates,
                            responses,
                            weights_to_use,
                            progress,
                        )
                    } else {
                        stochastic_dominance::uncensored::<Increasing>(
                            covariates,
                            responses,
                            weights_to_use,
                            progress,
                        )
                    }
                }
                StochasticOrder::HazardRateOrder => {
                    if decreasing {
                        hazard_rate_order::uncensored::<Decreasing>(
                            covariates,
                            responses,
                            weights_to_use,
                            settings.epsilon,
                            progress,
                        )
                    } else {
                        hazard_rate_order::uncensored::<Increasing>(
                            covariates,
                            responses,
                            weights_to_use,
                            settings.epsilon,
                            progress,
                        )
                    }
                }
            },
        }?;

        let AlgorithmOutput {
            mut cdfs,
            unique_covariates,
            thresholds,
        } = output;
        transpose(&mut cdfs, thresholds.len(), unique_covariates.len());

        Ok(Self {
            increasing: !decreasing,
            cdfs,
            covariates: unique_covariates,
            thresholds,
            quality_indicators: QualityIndicators {
                epsilon: settings.epsilon,
            },
        })
    }

    fn interpolate_covariate<'a>(
        &'a self,
        covariate: Self::Covariate<'_>,
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

    fn thresholds(&self) -> &[f64] {
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

pub struct GridPredictor<'a, I> {
    covariate_iter: I,
    state: GridPredictorState<'a>,
    thresholds: Vec<ResponseCoordinate>,
    threshold_index: usize,
}

impl<I: Iterator<Item = f64>> Iterator for GridPredictor<'_, I> {
    type Item = f64;

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

impl<I: ExactSizeIterator<Item = f64>> ExactSizeIterator for GridPredictor<'_, I> {}

impl Fit {
    pub fn predict_grid<I: IntoIterator<Item = f64>>(
        &self,
        covariates: I,
        thresholds: impl IntoIterator<Item = f64>,
    ) -> GridPredictor<'_, I::IntoIter> {
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

    /// An empty instance represents a sub-CDF that stays at value 0.0, always.
    fn is_empty(&self) -> bool {
        let result = self.cdfs.is_empty();
        debug_assert_eq!(self.covariates.is_empty(), result);
        debug_assert_eq!(self.thresholds.is_empty(), result);
        result
    }
}

#[derive(Debug, PartialEq)]
pub struct AlgorithmOutput {
    pub cdfs: Vec<f64>,
    pub unique_covariates: Vec<f64>,
    pub thresholds: Vec<f64>,
}

#[derive(Clone)]
pub struct Partition<W, V> {
    /// Right boundary of the partition, excluding. So partition is [previous, ..., last], index
    pub index: usize,
    pub weight: W,
    pub value: V,
}

pub type WeightedPartition = Partition<f64, f64>;

/// Information that is aggregated per covariate, like the total weight of all observations for a
/// covariate.
#[derive(Clone, Debug, PartialEq)]
pub struct CovariateStatistic {
    /// Total weight of this covariate, must be positive.
    pub weight: f64,
    /// Cumulative weight of this covariate and all that came before it (last covariate statistic
    /// then has the total weight of all observations).
    pub cumulative_weight: f64,
}

/// Input required to run the algorithm.
#[derive(Clone, Debug, PartialEq)]
pub struct AlgorithmContext<S> {
    /// Holds `n` items sorted by response from low to high, uncensored < censored, and covariate
    pub observations: Vec<Observation<usize, f64, S>>,
    /// Holds information about each covariate
    pub covariate_statistics: Vec<CovariateStatistic>,
    /// Thresholds
    pub unique_responses: Vec<f64>,
    /// Covariates
    pub unique_covariates: Vec<f64>,
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
    pub epsilon: f64,
}

#[cfg(test)]
mod test {
    use crate::IsotonicDistributionalRegressionFit;
    use crate::total_order::{Fit, QualityIndicators};

    #[test]
    fn test_predict() {
        let fit = Fit {
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
}
