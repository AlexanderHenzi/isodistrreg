pub mod functionals;
pub mod hazard_rate_order;
mod prediction;
mod preprocessing;
mod routines;
pub mod stochastic_dominance;
mod structures;
mod tonic_regression;

pub use functionals::algorithm;
pub use prediction::{CovariateSearch, GridPredictorState, Interpolation};
pub use structures::{Config, Fit, QualityIndicators};
pub use tonic_regression::algorithm as tonic_regression;
pub use tonic_regression::algorithm_pre_sorted_deduplicated as tonic_regression_pre_sorted;
