use crate::error::Error;
use crate::preprocessing::validate;
use crate::structures::{Direction, ExecutionMode};
use itertools::{izip, multiunzip};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cmp::Ordering;

/// Assume all unique, no duplicates, equal weights
pub fn algorithm(
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weight: &[f64],
) -> Result<((usize, usize), Vec<f64>), Error> {
    validate(x.chunks_exact(1), y, Some(observed), Some(weight))?;

    // Check that covariate and responses are unique
    assert_eq!(
        {
            let mut clone = y.to_vec();
            clone.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            clone.dedup();
            clone.len()
        },
        y.len(),
    );
    assert_eq!(
        {
            let mut clone = x.to_vec();
            clone.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            clone.dedup();
            clone.len()
        },
        x.len(),
    );

    let AlgorithmInput {
        y: response,
        x: covariate,
        observed: censoring,
        n_covariate,
        n_threshold,
        ..
    } = format_censored(x, y, observed, weight)?;

    let n = y.len();

    // Contains: (x index, y index, censoring), sorted by y index
    let data: Vec<_> = izip!(covariate, response, censoring).collect();
    let compute_survival_r = |r: usize| {
        let mut data_sub_filter = data
            .iter()
            .copied()
            .filter(|&(i, _, _)| i >= r)
            .collect::<Vec<_>>();
        debug_assert!(data_sub_filter.iter().all(|&(i, _, _)| i >= r));

        let mut survival_sy = Vec::with_capacity(n - r);
        for s in (r..n).rev() {
            debug_assert!(data_sub_filter.iter().all(|&(i, _, _)| i <= s));
            let n_rs = s - r + 1;

            let survival_y = data_sub_filter
                .iter()
                .enumerate()
                .scan(1.0, |s, (i_rs, &(_, j, c))| {
                    if !c {
                        *s *= (n_rs - i_rs - 1) as f64 / (n_rs - i_rs) as f64;
                    }
                    Some((j, *s))
                })
                .collect::<Vec<_>>();
            survival_sy.push(survival_y);
            data_sub_filter.retain(|&(i, _, _)| i != s);
        }
        survival_sy
    };

    #[cfg(feature = "parallel")]
    let mut survival_rsy = (0..n)
        .into_par_iter()
        .map(compute_survival_r)
        .collect::<Vec<_>>();
    #[cfg(not(feature = "parallel"))]
    let mut survival_rsy = (0..n).map(compute_survival_r).collect::<Vec<_>>();
    assert!(
        survival_rsy
            .iter()
            .flatten()
            .flatten()
            .all(|(_, s)| (0.0..=1.0).contains(s))
    );

    // let by_yrs: Vec<Vec<Vec<_>>> = (0..n).map(|j|
    //     (0..n).map(|r|
    //         (r..n).map(|s| {
    //             let y_values = &S_rsy[r][n - 1 - s];
    //             let result = y_values.binary_search_by_key(&j, |&(y, _)| y);
    //             match result {
    //                 Ok(index) => y_values[index].1,
    //                 Err(0) => 1.0,
    //                 Err(index) => y_values[index - 1].1,
    //             }
    //         }).collect()
    //     ).collect()
    // ).collect();

    // function of y, i and r
    let mut inner_term = (0..n)
        .rev()
        .map(|y_index| {
            (0..n)
                .map(|i| {
                    let map_fn = |r: &mut Vec<Vec<(usize, f64)>>| {
                        (i..n)
                            .map(|s| {
                                let relevant = &mut r[n - 1 - s];
                                relevant.pop_if(|(yy_index, _)| *yy_index > y_index);

                                match relevant.last() {
                                    Some(&(yy_index, value)) => {
                                        debug_assert!(yy_index <= y_index);
                                        value
                                    }
                                    None => 1.0,
                                }
                            })
                            .fold(f64::INFINITY, f64::min)
                    };
                    #[cfg(feature = "parallel")]
                    {
                        survival_rsy[0..=i]
                            .par_iter_mut()
                            .map(map_fn)
                            .collect::<Vec<_>>()
                    }
                    #[cfg(not(feature = "parallel"))]
                    {
                        survival_rsy[0..=i]
                            .iter_mut()
                            .map(map_fn)
                            .collect::<Vec<_>>()
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    inner_term.reverse();
    assert!(inner_term.iter().all(|x| (0..n).all(|r| {
        (r..n)
            .map(|i| x[i][r])
            .collect::<Vec<_>>()
            .windows(2)
            .all(|w| w[0] <= w[1])
    })));
    assert!(inner_term.iter().flatten().flatten().all(|v| v.is_finite()));

    // function of y and i
    let outer_term = inner_term.into_iter().map(|matrix| {
        matrix
            .into_iter()
            .enumerate()
            .map(|(i, rs)| rs.into_iter().take(i + 1).fold(f64::NEG_INFINITY, f64::max))
    });
    assert!(outer_term.clone().flatten().all(|v| v.is_finite()));

    let cdf = outer_term.flat_map(|s| s.map(|x| 1.0 - x)).collect();

    Ok(((n_covariate, n_threshold), cdf))
}

/// Input required to run the algorithm.
///
/// The covariate is always sorted itself.
#[derive(Clone)]
struct AlgorithmInput {
    /// Size `n_threshold`, potentially duplicate indices in `0..n_threshold`, sorted by response.
    y: Vec<usize>,
    /// Size n, indices in `0..n_covariate` sorted by response.
    x: Vec<usize>,
    /// Size n, sorted by covariate.
    observed: Vec<bool>,

    /// Number of unique covariate values.
    n_covariate: usize,
    /// Number of unique response (threshold) values.
    n_threshold: usize,
}

fn format_censored(
    covariate: &[f64],
    response: &[f64],
    observed: &[bool],
    weights: &[f64],
) -> Result<AlgorithmInput, Error> {
    let n = validate(
        covariate.chunks_exact(1),
        response,
        Some(observed),
        Some(weights),
    )?;

    let mut x_sorted_xywo = izip!(covariate, response, weights, observed)
        .map(|(&x, &y, &w, &o)| (x, y, w, o))
        .collect::<Vec<_>>();
    x_sorted_xywo.sort_unstable_by(|l, r| {
        l.0.partial_cmp(&r.0)
            .unwrap()
            .then(l.3.cmp(&r.3).reverse())
            .then(l.1.total_cmp(&r.1))
    });

    let mut x_index = 0;
    let mut weights_grouped = vec![x_sorted_xywo[x_index].2];
    let mut x = vec![x_index];
    let mut from_covariate_index_to_data_index = vec![0];
    for i in 1..n {
        if x_sorted_xywo[i - 1].0 == x_sorted_xywo[i].0 {
            weights_grouped[x_index] += x_sorted_xywo[i].2;
        } else {
            x_index += 1;
            weights_grouped.push(x_sorted_xywo[i].2);
            from_covariate_index_to_data_index.push(i);
        }
        x.push(x_index);
    }
    debug_assert_eq!(x.len(), n);
    let n_covariate = weights_grouped.len();
    debug_assert_eq!(from_covariate_index_to_data_index.len(), n_covariate);

    let mut y_sorted_xywo = izip!(x_sorted_xywo, x)
        .map(|((_, y, _w, o), x_index)| (x_index, y, o))
        .collect::<Vec<_>>();
    y_sorted_xywo.sort_by(|l, r| l.1.partial_cmp(&r.1).unwrap());
    let (response_groups, response_float, censoring): (Vec<_>, Vec<f64>, _) =
        multiunzip(y_sorted_xywo);

    let mut y_index = 0;
    let mut response = vec![y_index];
    for j in 1..n {
        if response_float[j - 1] < response_float[j] {
            y_index += 1;
        }
        response.push(y_index);
    }
    let n_threshold = y_index + 1;
    debug_assert_eq!(response.len(), n);

    Ok(AlgorithmInput {
        y: response,
        x: response_groups,
        observed: censoring,
        n_covariate,
        n_threshold,
    })
}

pub fn algorithm_single<D: Direction, EM: ExecutionMode>(
    threshold: f64,
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weights: &[f64],
) -> Result<Vec<f64>, Error> {
    validate(x.chunks_exact(1), y, Some(observed), Some(weights))?;

    let PreprocessingResult {
        below_threshold,
        from_covariate_index_to_data_index,
        cumulative_weight,
    } = format_censored_single(threshold, x, y, observed, weights)?;

    let n = from_covariate_index_to_data_index.len();
    assert_eq!(cumulative_weight.len(), n);
    assert!(from_covariate_index_to_data_index.iter().all(|&i| i < n));

    let compute_survival_row = |r: usize| {
        (r..n)
            .map(|s| {
                let start_index = from_covariate_index_to_data_index[r];
                let end_index = if s < n - 1 {
                    from_covariate_index_to_data_index[s + 1]
                } else {
                    below_threshold.len()
                };

                let mut data_subset = Vec::from(&below_threshold[start_index..end_index]);
                data_subset.sort_unstable_by(|&(y_l, o_l, _), &(y_r, o_r, _)| {
                    y_l.total_cmp(&y_r).then(o_l.cmp(&o_r).reverse())
                });
                let mut total_weight = if r > 0 {
                    cumulative_weight[s] - cumulative_weight[r - 1]
                } else {
                    cumulative_weight[s]
                };

                data_subset.iter().fold(1.0, |mut s, &(_, o, w)| {
                    if o {
                        s *= 1.0 - w / total_weight;
                    }
                    total_weight -= w;
                    s
                })
            })
            .collect()
    };

    #[cfg(feature = "parallel")]
    let survival_rs: Vec<Vec<_>> = if EM::PARALLEL {
        (0..n).into_par_iter().map(compute_survival_row).collect()
    } else {
        (0..n).map(compute_survival_row).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let survival_rs: Vec<Vec<_>> = (0..n).map(compute_survival_row).collect();
    assert!(
        survival_rs
            .iter()
            .flatten()
            .all(|s| (0.0..=1.0).contains(s))
    );

    type MinOrMax = fn(f64, f64) -> f64;
    let (outer_direction, inner_direction): (MinOrMax, MinOrMax) = match D::FORBIDDEN_ORDERING {
        Ordering::Less => (f64::max, f64::min),
        Ordering::Equal => panic!(),
        Ordering::Greater => (f64::min, f64::max),
    };

    let extremes = (0..n)
        .map(|i| {
            (0..=i)
                .map(|r| {
                    (i..n)
                        .map(|s| survival_rs[r][s - r])
                        .reduce(inner_direction)
                        .unwrap()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let results = extremes
        .into_iter()
        .map(|r_values| r_values.into_iter().reduce(outer_direction).unwrap())
        .map(|s| 1.0 - s)
        .collect::<Vec<_>>();

    Ok(results)
}

struct PreprocessingResult {
    below_threshold: Vec<(f64, bool, f64)>,
    from_covariate_index_to_data_index: Vec<usize>,
    cumulative_weight: Vec<f64>,
}

fn format_censored_single(
    threshold: f64,
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    weights: &[f64],
) -> Result<PreprocessingResult, Error> {
    let n = validate(x.chunks_exact(1), y, Some(observed), Some(weights))?;

    let mut x_sorted_xycw = izip!(x, y, observed, weights)
        .map(|(&x, &y, &c, &w)| (x, y, c, w))
        .collect::<Vec<_>>();
    x_sorted_xycw.sort_unstable_by(|&(x_l, _, _, _), &(x_r, _, _, _)| x_l.total_cmp(&x_r));

    let mut below_threshold = vec![];

    let first = x_sorted_xycw[0];
    let mut weights_grouped = vec![first.3];
    let mut from_covariate_index_to_data_index = vec![below_threshold.len()];
    if x_sorted_xycw[0].1 <= threshold {
        below_threshold.push((first.1, first.2, first.3));
    }
    for i in 1..n {
        let current = x_sorted_xycw[i];

        if x_sorted_xycw[i - 1].0 == current.0 {
            // Same covariate
            weights_grouped[from_covariate_index_to_data_index.len() - 1] += current.3;
        } else {
            // New covariate
            below_threshold[*from_covariate_index_to_data_index.last().unwrap()..]
                .sort_unstable_by(|&(y_l, c_l, _), &(y_r, c_r, _)| {
                    y_l.total_cmp(&y_r).then(c_l.cmp(&c_r))
                });
            weights_grouped.push(current.3);
            from_covariate_index_to_data_index.push(below_threshold.len());
        }

        if current.1 <= threshold {
            below_threshold.push((current.1, current.2, current.3));
        }
    }

    let cumulative_weight: Vec<_> = weights_grouped
        .into_iter()
        .scan(0.0, |total, value| {
            *total += value;
            Some(*total)
        })
        .collect();

    debug_assert_eq!(
        from_covariate_index_to_data_index.len(),
        cumulative_weight.len()
    );

    Ok(PreprocessingResult {
        below_threshold,
        from_covariate_index_to_data_index,
        cumulative_weight,
    })
}

#[cfg(test)]
mod test {
    use crate::structures::{Decreasing, Serial};
    use crate::total_order::stochastic_dominance::censored::nonrecursive::algorithm_single;

    #[test]
    fn test_single() {
        assert_eq!(
            algorithm_single::<Decreasing, Serial>(
                3.5,
                &[1.0, 2.0, 3.0, 4.0],
                &[2.0, 1.0, 4.0, 3.0],
                &[false, true, true, true],
                &[1.0; 4],
            )
            .unwrap(),
            vec![5.0 / 8.0, 5.0 / 8.0, 0.5, 0.5],
        );
    }
}
