use crate::error::Error;
use crate::progress::ProgressTracker;
use crate::routines::empirical_cdf;
use crate::structures::Direction;
use crate::total_order::preprocessing::preprocess_uncensored;
use crate::total_order::stochastic_dominance::routines;
use crate::total_order::structures::{AlgorithmContext, AlgorithmOutput};
use std::iter::repeat_n;

pub fn algorithm<D: Direction>(
    x: &[f64],
    y: &[f64],
    weight: &[f64],
    progress: &dyn ProgressTracker,
) -> Result<AlgorithmOutput, Error> {
    let AlgorithmContext {
        observations,
        covariate_statistics,
        unique_responses,
        unique_covariates,
    } = preprocess_uncensored(x, y, weight);
    progress.set_total(unique_responses.len());

    if unique_responses.len() == 1 {
        // Single threshold -> all one's
        return Ok(AlgorithmOutput {
            cdfs: vec![1.0; unique_covariates.len()],
            unique_covariates,
            thresholds: unique_responses,
        });
    }
    if covariate_statistics.len() == 1 {
        // Single covariate -> a single empirical cdf
        // Observations have been deduplicated so responses are unique
        return Ok(AlgorithmOutput {
            cdfs: empirical_cdf(observations.iter().copied(), covariate_statistics[0].weight),
            unique_covariates,
            thresholds: unique_responses,
        });
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
        &observations,
        &covariate_statistics,
        |()| false,
        &mut partitions_to_store,
        &mut cdfs,
        progress,
    );

    // Final threshold
    cdfs.extend(repeat_n(1.0, unique_covariates.len()));
    progress.increment();

    cdfs.shrink_to_fit();

    Ok(AlgorithmOutput {
        cdfs,
        unique_covariates,
        thresholds: unique_responses,
    })
}
