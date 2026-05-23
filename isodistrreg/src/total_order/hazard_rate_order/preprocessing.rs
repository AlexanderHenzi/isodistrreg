use crate::Float;
use crate::total_order::preprocessing;
use crate::total_order::structures::AlgorithmContext;

pub fn preprocess_censored<X: Float, Y: Float, W: Float>(
    x: &[X],
    y: &[Y],
    observed: &[bool],
    weight: &[W],
) -> AlgorithmContext<X, Y, bool> {
    let mut context = preprocessing::preprocess(x, y, |i| observed[i], weight);

    context.unique_responses = context
        .observations
        .first()
        .into_iter()
        .chain(
            context
                .observations
                .array_windows()
                .filter_map(|[left, right]| (right.y != left.y).then_some(right)),
        )
        .filter(|o| o.observed)
        .map(|o| o.y)
        .collect();

    context
}
