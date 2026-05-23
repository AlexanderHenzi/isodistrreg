//! Shared noise-clamping tolerance for the f32 post-preprocessing algorithms.
//!
//! Several total-order algorithms (the SD-censored fast path and both hazard-rate paths)
//! run the entire inner loop in f32 and need to distinguish "this weight is exactly zero
//! up to f32 round-off" from "this weight is genuinely zero". The conservative bound
//! below works for all of them; per-call-site interpretation lives next to each call.

/// Dynamic f32-noise floor for sums of `n_total` weighted observations.
///
/// The algorithms maintain running sums and differences of f32 weights. A Wilkinson-style
/// worst-case error bound for a length-`n` sum (or for a cum-diff plus a running
/// subtraction) is `3·n·u_32` where `u_32 = 2^-24` is the formal unit roundoff. We
/// double it for the "sum minus a running subtraction" pattern (≈ `6·n·u_32`) and pick
/// 8× the unit roundoff for headroom:
///
///   `epsilon = 8 · n · u_32 ≈ n · 10⁻⁶`
///
/// In code we multiply by `f32::EPSILON` (= `2^-23` = `2 · u_32`), so the literal here is
/// `8 · n · f32::EPSILON`; the 8× margin absorbs the factor-of-2 between `f32::EPSILON`
/// and `u_32` either way.
///
/// At `n ≤ 10⁵` this stays well below 1 (the weight of a single unit-weight observation),
/// so no real cell is mistakenly skipped. At `n ≥ 10⁷` in f32 the algorithm has lost
/// useful precision regardless; this threshold clamps it to keep it numerically alive,
/// but the output is no longer accurate at any threshold.
#[inline]
pub fn weight_noise_floor(n_total: usize) -> f32 {
    // `f32::EPSILON` is the gap from 1.0 to the next f32 (2^-23), i.e. 2× the formal
    // numerical-analysis "unit roundoff" u_32 = 2^-24. The 8× margin above absorbs that
    // factor-of-2 either way.
    8.0 * (n_total as f32) * f32::EPSILON
}
