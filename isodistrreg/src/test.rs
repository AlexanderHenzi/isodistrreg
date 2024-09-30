use approx::relative_eq;

pub fn is_relative_eq_vec(result: &[f64], expected: &[f64]) -> bool {
    if result.len() != expected.len() {
        false
    } else {
        result
            .iter()
            .zip(expected.iter())
            .all(|(l, r)| relative_eq!(l, r))
    }
}
