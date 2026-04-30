#[cfg(test)]
mod definition;
mod fast;
mod nonrecursive;
mod preprocessing;
mod propagate_bounds;
mod structures;

#[cfg(test)]
pub(crate) use definition::algorithm as algorithm_definition;
pub use fast::algorithm;
pub use nonrecursive::algorithm as algorithm_nonrecursive;
pub use nonrecursive::algorithm_single as algorithm_nonrecursive_single;
#[cfg(test)]
pub(crate) use preprocessing::preprocess;

#[cfg(test)]
mod test;
