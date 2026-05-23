mod censored;
mod preprocessing;
mod routines;
mod uncensored;

pub use censored::algorithm as censored;
pub use preprocessing::preprocess_censored;
pub use uncensored::algorithm as uncensored;

// Re-export the shared uncensored preprocessor used by the hazard-rate uncensored algorithm so
// callers can build an `AlgorithmContext` to pass into `uncensored`.
pub use crate::total_order::preprocessing::preprocess_uncensored;
