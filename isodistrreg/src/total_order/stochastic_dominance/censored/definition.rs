use crate::error::Error;
use crate::structures::Direction;
use crate::total_order::stochastic_dominance::censored::preprocess;
use crate::total_order::stochastic_dominance::censored::preprocessing::PreProcessingResult;
use crate::total_order::structures::AlgorithmOutput;

pub(crate) fn algorithm<D: Direction>(
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weight: &[f64],
) -> Result<AlgorithmOutput, Error> {
    let PreProcessingResult {
        context,
        unique_covariates,
        thresholds,
    } = preprocess(x, y, observed, weight)?;

    let cdfs: Vec<_> = (0..thresholds.len())
        .map(|threshold| {
            let mut survivals = (0..context.n_covariate())
                .map(|r| {
                    (r..context.n_covariate())
                        .map(|s| {
                            let data = context
                                .observations
                                .iter()
                                .filter(|o| r <= o.x && o.x <= s)
                                .collect::<Vec<_>>();
                            let mut survival = 1.0;
                            for (i_rs, observation) in data.iter().enumerate() {
                                if threshold.lt(&observation.y) {
                                    break;
                                }

                                if observation.observed {
                                    survival *= 1.0
                                        - observation.weight
                                            / data[i_rs..].iter().map(|o| o.weight).sum::<f64>();
                                }
                            }

                            survival
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();

            // Clip the estimators
            fn clip(value: f64, left: f64, right: f64) -> f64 {
                let minimum = left.min(right);
                let maximum = left.max(right);

                value.min(maximum).max(minimum)
            }
            for r in 0..context.n_covariate() - 1 {
                for s in r + 1..context.n_covariate() {
                    for k in r..s {
                        // All bounds are inclusive
                        let left_start = r;
                        let left_end = k;
                        let right_start = k + 1;
                        let right_end = s;

                        survivals[r][s - r] = clip(
                            survivals[r][s - r],
                            survivals[left_start][left_end - left_start],
                            survivals[right_start][right_end - right_start],
                        )
                    }
                }
            }

            // Each threshold is an isotonic regression in the opposite direction of the ordering of
            // the survival quantities. Normal IDR (increasing) is for a single threshold min-max
            // (decreasing), but because we work with survival quantities here, it's max-min instead
            // for S-IDR increasing.
            type Red = fn(f64, f64) -> f64;
            let (outer_reduction, inner_reduction): (Red, Red) = match D::IS_INCREASING {
                true => (f64::max, f64::min),
                false => (f64::min, f64::max),
            };

            // Take inner extreme (min for increasing)
            for r in 0..context.n_covariate() {
                for s in (r..context.n_covariate() - 1).rev() {
                    let s_data = s - r;

                    survivals[r][s_data] =
                        inner_reduction(survivals[r][s_data], survivals[r][s_data + 1]);
                }
            }
            // Take outer extreme (max for decreasing)
            let mut maximums = Vec::with_capacity(context.n_covariate());
            for i in 0..context.n_covariate() {
                // Starting value for the max / min
                let mut candidate = match D::IS_INCREASING {
                    true => f64::NEG_INFINITY,
                    false => f64::INFINITY,
                };
                for r in 0..=i {
                    let i_data = i - r;
                    candidate = outer_reduction(candidate, survivals[r][i_data]);
                }
                maximums.push(candidate);
            }

            maximums
        })
        .flatten()
        .map(|s| 1.0 - s)
        .collect();

    Ok(AlgorithmOutput {
        cdfs,
        unique_covariates,
        thresholds,
    })
}
