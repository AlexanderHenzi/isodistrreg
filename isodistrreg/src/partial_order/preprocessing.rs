use crate::error::Error;
use crate::partial_order::ExtendedAlgorithmContext;
use crate::partial_order::structures::{AlgorithmContext, CovariateGroups, PartialOrder};
use crate::routines::{argsort_unstable_by, lexicographic_cmp};
use crate::structures::Increasing;
use itertools::Itertools;
use std::cmp::Ordering;
use std::collections::HashMap;

pub fn validate(
    covariates: &[f64],
    responses: &[f64],
    censoring: Option<&[bool]>,
    weights: Option<&[f64]>,
    covariate_order: &CovariateGroups,
) -> Result<usize, Error> {
    if !covariate_order.is_consistent() {
        return Err(Error::CovariateOrderInconsistency);
    }

    let n = crate::preprocessing::validate(
        covariates.chunks_exact(covariate_order.dimension),
        responses,
        censoring,
        weights,
    )?;

    if covariates.len() != n * covariate_order.dimension {
        return Err(Error::IncompatibleShapes {
            covariate_len: covariates.len(),
            response_len: responses.len(),
            weight_len: weights.map(|slice| slice.len()),
            y_observed_len: censoring.map(|slice| slice.len()),
        });
    }

    Ok(n)
}

#[must_use]
pub fn preprocess_uncensored(
    covariates: &[f64],
    responses: &[f64],
    weights: &[f64],
    covariate_groups: &CovariateGroups,
) -> AlgorithmContext {
    // Convert the covariates such that they become comparable component-wise
    let covariates: Vec<_> = covariates
        .chunks_exact(covariate_groups.dimension)
        .flat_map(|covariate| unify_group_orders(covariate, covariate_groups))
        .collect();
    // Aggregate duplicates: unique covariates, and sorted (covariate, response, weight) tuples
    group_by_covariate_uncensored(&covariates, responses, weights, covariate_groups.dimension)
}

#[must_use]
pub fn preprocess_censored(
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weight: &[f64],
    covariate_groups: &CovariateGroups,
) -> ExtendedAlgorithmContext {
    // Convert the covariates such that they become comparable component-wise
    let x_converted: Vec<_> = x
        .chunks_exact(covariate_groups.dimension)
        .flat_map(|covariate| unify_group_orders(covariate, covariate_groups))
        .collect();
    // Aggregate duplicates: unique covariates, and sorted (covariate, response, weight) tuples
    group_by_covariate_censored(
        &x_converted,
        y,
        observed,
        weight,
        covariate_groups.dimension,
    )
}

/// Modify the covariate data such that the covariates can be compared component-wise.
///
/// Groups of covariates should be treated together if their values are interchangeable. We simplify
/// downstream logic by converting the covariate groups such that their ordering can be determined
/// by a component-wise comparison.
pub fn unify_group_orders(
    covariate: &[f64],
    covariate_groups: &CovariateGroups,
) -> impl Iterator<Item = f64> {
    let mut result = Vec::with_capacity(covariate.len());

    // Re-organize (within each group) the values such that we can compare them component-wise
    // after, independent of the initial order constraint on the group.
    for (order, columns) in &covariate_groups.groups {
        let group_elements = columns.iter().map(|&c| covariate[c]);
        match order {
            PartialOrder::ComponentWise => {
                // Do not modify the group elements - we compare them element-wise as they are
                result.extend(group_elements);
            }
            PartialOrder::StochasticDominance => {
                // Sort the elements within the group descending
                let previous_len = result.len();
                result.extend(group_elements);
                result[previous_len..].sort_unstable_by(|l, r| r.total_cmp(l));
            }
            PartialOrder::IncreasingConvexOrder => {
                // Sort the elements within the group descending
                let previous_len = result.len();
                result.extend(group_elements);
                result[previous_len..].sort_unstable_by(|l, r| r.total_cmp(l));
                // ... and compute their cumulative sum
                for k in previous_len + 1..result.len() {
                    result[k] += result[k - 1];
                }
            }
        }
    }

    result
        .into_iter()
        .chain(covariate_groups.ungrouped.iter().map(|&c| covariate[c]))
}

/// Group rows of X; for each unique covariate row keep all responses and sum of weights.
fn group_by_covariate_uncensored(
    covariates: &[f64],
    responses: &[f64],
    weights: &[f64],
    covariate_dimension: usize,
) -> AlgorithmContext {
    let n_max = responses.len();
    debug_assert_eq!(covariates.len(), n_max * covariate_dimension);

    // Sort by covariate, response and censoring to copy deduplicating, building a covariate index

    let get_cdf = |i| &covariates[i * covariate_dimension..(i + 1) * covariate_dimension];
    let covariate_order = argsort_unstable_by::<Increasing, _>(
        |a, b| {
            lexicographic_cmp(get_cdf(a), get_cdf(b)).then(responses[a].total_cmp(&responses[b]))
        },
        n_max,
    );

    let n_initial_covariates = covariates.len() / covariate_dimension;
    // Unique, sorted covariates
    let mut sorted_covariates = Vec::with_capacity(n_initial_covariates);
    // Observation weight aggregated by covariate
    let mut total_weight_per_covariate = Vec::with_capacity(n_initial_covariates);
    // Response value sorted by covariate (first) and response (second), their combination is unique
    let mut sorted_responses = Vec::with_capacity(covariates.len());
    // Weight values aggregated by sorted covariate (first) and response (second)
    let mut sorted_weights = Vec::with_capacity(covariates.len());
    // For each observation, the index in which the covariate is stored, by the same order
    let mut covariate_indices = Vec::with_capacity(n_max);

    // First item
    let read_index = covariate_order[0];
    sorted_covariates.extend_from_slice(get_cdf(read_index));
    total_weight_per_covariate.push(0.0);
    sorted_responses.push(responses[read_index]);
    sorted_weights.push(weights[read_index]);
    covariate_indices.push(0);

    // Remaining items
    for read_index in covariate_order.into_iter().skip(1) {
        let nr_unique = sorted_covariates.len() / covariate_dimension;
        let last = &sorted_covariates[(nr_unique - 1) * covariate_dimension..];
        if get_cdf(read_index) != last {
            // a new covariate

            // from last iteration, collect weight
            *total_weight_per_covariate.last_mut().unwrap() += sorted_weights.last().unwrap();

            // this iteration, add the new data
            sorted_covariates.extend_from_slice(get_cdf(read_index));
            total_weight_per_covariate.push(0.0);
            sorted_responses.push(responses[read_index]);
            sorted_weights.push(weights[read_index]);
            covariate_indices.push(nr_unique);
        } else if responses[read_index] != *sorted_responses.last().unwrap() {
            // same covariate, but new response

            // from last iteration, collect weight
            *total_weight_per_covariate.last_mut().unwrap() += sorted_weights.last().unwrap();

            // this iteration, add the new data
            sorted_responses.push(responses[read_index]);
            sorted_weights.push(weights[read_index]);
            covariate_indices.push(nr_unique - 1); // still the old index
        } else {
            // same covariate and response, just increment the weight
            *sorted_weights.last_mut().unwrap() += weights[read_index];
        }

        let nr_covariates = total_weight_per_covariate.len();
        assert_eq!(sorted_covariates.len(), nr_covariates * covariate_dimension);
        assert_eq!(total_weight_per_covariate.len(), nr_covariates);
        let nr_observations = sorted_responses.len();
        assert_eq!(sorted_responses.len(), nr_observations);
        assert_eq!(sorted_weights.len(), nr_observations);
        assert_eq!(covariate_indices.len(), nr_observations);
    }
    *total_weight_per_covariate.last_mut().unwrap() += sorted_weights.last().unwrap();
    debug_assert!(
        (total_weight_per_covariate.iter().sum::<f64>() - sorted_weights.iter().sum::<f64>()).abs()
            < 1e-6
    );
    // Gets stored in the final Fit struct, so let's not waste space
    sorted_covariates.shrink_to_fit();

    for (c, tw) in total_weight_per_covariate.iter().enumerate() {
        let check = covariate_indices
            .iter()
            .zip(sorted_weights.iter())
            .filter(|&(&cc, _)| cc == c)
            .map(|(_, &w)| w)
            .sum::<f64>();
        debug_assert!((check - tw).abs() < 1e-6);
    }

    let order = argsort_unstable_by::<Increasing, _>(
        |l, r| sorted_responses[l].total_cmp(&sorted_responses[r]),
        sorted_responses.len(),
    );
    let (responses, weights, covariate_indices) = order
        .into_iter()
        .map(|i| (sorted_responses[i], sorted_weights[i], covariate_indices[i]))
        .collect::<(Vec<_>, Vec<_>, Vec<_>)>();
    let thresholds = responses.iter().copied().dedup().collect();

    for (c, tw) in total_weight_per_covariate.iter().enumerate() {
        let check = covariate_indices
            .iter()
            .zip(weights.iter())
            .filter(|&(&cc, _)| cc == c)
            .map(|(_, &w)| w)
            .sum::<f64>();
        debug_assert!((check - tw).abs() < 1e-6);
    }

    // Gets stored in the final Fit struct, so let's not waste space
    sorted_covariates.shrink_to_fit();

    AlgorithmContext {
        x: sorted_covariates,
        x_weight: total_weight_per_covariate,
        y: responses,
        weight: weights,
        x_indices: covariate_indices,
        thresholds,
    }
}

#[derive(Clone, Copy, Debug)]
struct ObsTmp {
    /// row in the *original* covariate matrix
    x_row: usize,
    /// index into `thresholds` (a.k.a. \( \tau \))
    threshold: usize,
    observed: bool,
    weight: f64,
}

/// 1) sorts by (response, censoring, covariate), discards initial censored rows,
///    builds unique response thresholds and stores threshold indices;
/// 2) deduplicates covariates (multi-dim) by value (lexicographic);
/// 3) re-sorts by (threshold/response, censoring) and deduplicates again.
fn group_by_covariate_censored(
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weight: &[f64],
    dimension: usize,
) -> ExtendedAlgorithmContext {
    let n = y.len();
    debug_assert_eq!(observed.len(), n);
    debug_assert_eq!(weight.len(), n);
    debug_assert_eq!(x.len(), n * dimension);

    let get_cdf = |i| &x[i * dimension..(i + 1) * dimension];

    // -------------------------------------------------------------------------
    // 1) Sort by (response, censoring, covariate). Drop all initial censored obs.
    //    Build thresholds and store threshold index per observation.
    // -------------------------------------------------------------------------

    let response_order = argsort_unstable_by::<Increasing, _>(
        |i, j| {
            y[i].total_cmp(&y[j])
                .then(observed[i].cmp(&observed[j]).reverse())
                .then_with(|| lexicographic_cmp(get_cdf(i), get_cdf(j)))
        },
        n,
    );

    // Find first uncensored in the (response, censoring, covariate)-sorted order.
    let mut first_uncensored_pos = None;
    for (pos, &idx) in response_order.iter().enumerate() {
        if observed[idx] {
            first_uncensored_pos = Some(pos);
            break;
        }
    }

    // If everything is censored, all thresholds are empty and we discard everything.
    let Some(start) = first_uncensored_pos else {
        return ExtendedAlgorithmContext {
            x: Vec::with_capacity(0),
            y: Vec::with_capacity(0),
            y_observed: Vec::with_capacity(0),
            weight: Vec::with_capacity(0),
            x_unique: Vec::with_capacity(0),
            x_weight: Vec::with_capacity(0),
            thresholds: Vec::with_capacity(0),
        };
    };

    let mut thresholds: Vec<f64> = Vec::with_capacity(n - start);
    let mut obs: Vec<ObsTmp> = Vec::with_capacity(n - start);

    // Helper: push (or merge) an observation in the current sorted stream.
    let push_or_merge = |o: ObsTmp, obs: &mut Vec<ObsTmp>| {
        if let Some(last) = obs.last_mut()
            && last.x_row == o.x_row
            && last.threshold == o.threshold
            && last.observed == o.observed
        {
            last.weight += o.weight;
        } else {
            obs.push(o);
        }
    };

    for &data_idx in &response_order[start..] {
        // Update thresholds only on *uncensored* response changes.
        if observed[data_idx] && thresholds.last().copied() != Some(y[data_idx]) {
            thresholds.push(y[data_idx]);
        }

        // Map censored obs to the neighboring (previous / lower-or-equal) uncensored threshold index.
        // This matches the 1D reference implementation: censored -> current last threshold.
        let threshold_index = thresholds.len().saturating_sub(1);

        // Note: after discarding initial censored observations, `thresholds` is non-empty here.
        debug_assert!(!thresholds.is_empty());

        push_or_merge(
            ObsTmp {
                x_row: data_idx,
                threshold: threshold_index,
                observed: observed[data_idx],
                weight: weight[data_idx],
            },
            &mut obs,
        );
    }

    // -------------------------------------------------------------------------
    // 2) Sort by covariate and deduplicate covariate rows by value.
    //    Build mapping original cov_row -> unique covariate index.
    // -------------------------------------------------------------------------

    // Only consider covariate rows that survived step 1.
    let mut cov_rows: Vec<usize> = obs.iter().map(|o| o.x_row).collect();
    cov_rows.sort_unstable();
    cov_rows.dedup();

    cov_rows.sort_unstable_by(|&i, &j| lexicographic_cmp(get_cdf(i), get_cdf(j)));

    let mut unique_covariates: Vec<f64> = Vec::with_capacity(cov_rows.len() * dimension);
    let mut cov_row_to_unique: HashMap<usize, usize> = HashMap::with_capacity(cov_rows.len());

    let mut last_row: Option<usize> = None;
    let mut current_unique = usize::MAX;

    for &row in &cov_rows {
        if last_row.is_none()
            || lexicographic_cmp(get_cdf(row), get_cdf(last_row.unwrap())) != Ordering::Equal
        {
            unique_covariates.extend_from_slice(&x[row * dimension..(row + 1) * dimension]);
            current_unique = unique_covariates.len() / dimension - 1;
            last_row = Some(row);
        }
        cov_row_to_unique.insert(row, current_unique);
    }

    // Rewrite observations to refer to the unique covariate index (not the original row).
    for o in &mut obs {
        o.x_row = *cov_row_to_unique.get(&o.x_row).unwrap();
    }

    // -------------------------------------------------------------------------
    // 3) Final sort by (threshold/response, censoring) again, then deduplicate.
    //    Also compute total weight per unique covariate.
    // -------------------------------------------------------------------------

    let mut order: Vec<usize> = (0..obs.len()).collect();
    order.sort_unstable_by(|&l, &r| {
        obs[l]
            .threshold
            .cmp(&obs[r].threshold)
            .then(obs[l].observed.cmp(&obs[r].observed).reverse())
            .then(obs[l].x_row.cmp(&obs[r].x_row))
    });

    let mut out_threshold_idx: Vec<usize> = Vec::with_capacity(obs.len());
    let mut out_observed: Vec<bool> = Vec::with_capacity(obs.len());
    let mut out_weight: Vec<f64> = Vec::with_capacity(obs.len());
    let mut out_cov_idx: Vec<usize> = Vec::with_capacity(obs.len());

    for &k in &order {
        let o = obs[k];
        if let Some(last_i) = out_threshold_idx.len().checked_sub(1)
            && out_threshold_idx[last_i] == o.threshold
            && out_observed[last_i] == o.observed
            && out_cov_idx[last_i] == o.x_row
        {
            out_weight[last_i] += o.weight;
        } else {
            out_threshold_idx.push(o.threshold);
            out_observed.push(o.observed);
            out_weight.push(o.weight);
            out_cov_idx.push(o.x_row);
        }
    }

    let nr_unique_cov = unique_covariates.len() / dimension;
    let mut covariate_weight = vec![0.0; nr_unique_cov];
    for (&c, &w) in out_cov_idx.iter().zip(out_weight.iter()) {
        covariate_weight[c] += w;
    }

    unique_covariates.shrink_to_fit();

    ExtendedAlgorithmContext {
        x: out_cov_idx,
        y: out_threshold_idx,
        y_observed: out_observed,
        weight: out_weight,
        x_unique: unique_covariates,
        x_weight: covariate_weight,
        thresholds,
    }
}
