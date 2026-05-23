use approx::AbsDiffEq;
use approx::relative_eq;

/// Element-wise approximate-equality check that auto-adapts to f32 or f64 via the `AbsDiffEq`
/// default epsilon, so test bodies that pass f32 cdf slices (now the production output type)
/// don't need to specify tolerances explicitly.
pub fn is_relative_eq_vec<T>(result: &[T], expected: &[T]) -> bool
where
    T: AbsDiffEq<Epsilon = T> + Copy + approx::RelativeEq,
{
    if result.len() != expected.len() {
        false
    } else {
        result
            .iter()
            .zip(expected.iter())
            .all(|(l, r)| relative_eq!(l, r))
    }
}
