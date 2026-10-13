#[cfg(test)]
mod definition;
#[allow(unsafe_code)]
mod fast;
mod nonrecursive;
#[allow(unsafe_code)]
mod propagate_bounds;
#[allow(unsafe_code)]
mod structures;

#[cfg(test)]
pub(crate) use definition::algorithm as algorithm_definition;
pub use fast::algorithm;
pub use nonrecursive::algorithm as algorithm_nonrecursive;
pub use nonrecursive::algorithm_single as algorithm_nonrecursive_single;

#[cfg(test)]
mod test;

#[cfg(test)]
mod scale_test;
