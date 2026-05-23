use crate::functionals::TotalCmp;

/// Float scalar types supported as the covariate (`X`) or threshold (`Y`) precision.
///
/// Extends `num_traits::Float` with [`TotalCmp`], which `num_traits::Float` does not expose
/// but every algorithm in this crate relies on (it gives a deterministic ordering on f32/f64
/// values, including NaN). Implementations are sealed to `f32` and `f64`.
///
/// `Default` is also required: `f32`/`f64` satisfy it trivially, and downstream algorithms
/// (notably the censored partial-order solver via `KaplanMeier<Y>`) need it without having
/// to spell out the conjunction at every call site.
pub trait Float: num_traits::Float + Default + TotalCmp + Send + Sync + 'static {}

impl Float for f32 {}

impl Float for f64 {}
