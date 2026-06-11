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

    if context.unique_responses.is_empty() {
        // No events at all: there are no thresholds, so the fit is the fully-empty
        // fit (a sub-CDF that is 0 everywhere), like the SD-censored route. Clear
        // the covariate grid too — `Fit::is_empty` requires covariates, thresholds
        // and cdfs to be empty together.
        context.observations = Vec::with_capacity(0);
        context.covariate_statistics = Vec::with_capacity(0);
        context.unique_covariates = Vec::with_capacity(0);
    }

    context
}
