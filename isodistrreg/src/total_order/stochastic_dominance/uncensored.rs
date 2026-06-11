use crate::Float;
use crate::progress::ProgressTracker;
use crate::routines::empirical_cdf;
use crate::structures::Direction;
use crate::total_order::stochastic_dominance::routines;
use crate::total_order::structures::AlgorithmContext;
use std::iter::repeat_n;

pub fn algorithm<D: Direction, X: Float, Y: Float>(
    context: &AlgorithmContext<X, Y, ()>,
    progress: &dyn ProgressTracker,
) -> Vec<f32> {
    let AlgorithmContext {
        observations,
        covariate_statistics,
        unique_responses,
        unique_covariates,
    } = context;
    progress.set_total(unique_responses.len());

    if observations.is_empty() {
        // Preprocessing dropped every observation (all had zero weight): the empty
        // fit, with no thresholds and no covariates.
        return Vec::with_capacity(0);
    }
    if unique_responses.len() == 1 {
        // Single threshold -> all one's
        return vec![1.0; unique_covariates.len()];
    }
    if covariate_statistics.len() == 1 {
        // Single covariate -> a single empirical cdf
        // Observations have been deduplicated so responses are unique
        return empirical_cdf(observations.iter().copied(), covariate_statistics[0].weight);
    }

    // TODO: Try asserting all we know is true about the input to allow the compiler to make more
    //  assumptions

    // Collects final estimate
    let mut cdfs = Vec::with_capacity(unique_responses.len() * unique_covariates.len());

    // Tracks which observation we're treating, sorted by response and censoring (and covariate)
    let mut data_index = 0;

    let mut partitions_to_store = Vec::with_capacity(unique_covariates.len());

    routines::accelerated_pava::<_, _, _, D>(
        &mut data_index,
        observations,
        covariate_statistics,
        |()| false,
        &mut partitions_to_store,
        &mut cdfs,
        progress,
    );

    // Final threshold
    cdfs.extend(repeat_n(1.0, unique_covariates.len()));
    progress.increment();

    cdfs.shrink_to_fit();

    cdfs
}
