//! Shared noise-clamping tolerance for the f32 post-preprocessing algorithms.
//!
//! Several total-order algorithms (the SD-censored fast path and both hazard-rate paths)
//! run the entire inner loop in f32 and need to distinguish "this weight is exactly zero
//! up to f32 round-off" from "this weight is genuinely zero". The conservative bound
//! below works for all of them; per-call-site interpretation lives next to each call.

/// Dynamic f32-noise floor for running sums and differences of observation weights
/// totalling `total_weight`.
///
/// The algorithms maintain running sums and differences of f32 weights. A Wilkinson-style
/// worst-case error bound for such accumulations is proportional to `u_32 · Σ|w_i|`
/// where `u_32 = 2^-24` is the formal unit roundoff — the error scales with the
/// magnitude of the weights, not with their count alone. We pick 8× the unit roundoff
/// for headroom:
///
///   `epsilon = 8 · Σw · u_32`
///
/// For unit weights this equals the count-based form `8 · n · u_32 ≈ n · 10⁻⁶`,
/// staying well below 1 (the weight of a single observation) for `n ≤ 10⁵`. Unlike the
/// count-based form it is invariant under rescaling all weights by a positive constant,
/// matching the scale invariance of the Kaplan-Meier and least-squares estimators the
/// kernels compute: externally normalized weights (e.g. summing to 1) no longer trip
/// the exhausted-weight shortcuts on cells that still hold real observations.
///
/// In code we multiply by `f32::EPSILON` (= `2^-23` = `2 · u_32`); the 8× margin absorbs
/// the factor-of-2 between `f32::EPSILON` and `u_32` either way.
#[inline]
pub fn weight_noise_floor(total_weight: f32) -> f32 {
    8.0 * total_weight * f32::EPSILON
}
