mod censored;
mod routines;
mod uncensored;

pub use censored::algorithm as censored;
pub use censored::algorithm_nonrecursive as censored_nonrecursive;
pub use censored::algorithm_nonrecursive_single as censored_nonrecursive_single;
pub(crate) use censored::preprocess as preprocess_censored;
pub use uncensored::algorithm as uncensored;

#[cfg(test)]
mod test;
