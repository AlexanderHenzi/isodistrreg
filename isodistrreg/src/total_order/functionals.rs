//! TODO: Experimental module with isotonic regressions for some functionals that violate the Cauchy
//!  mean-value property.

use crate::functionals::{CauchyMeanValueFunctional, Functional};
use crate::routines::median;
use crate::structures::Direction;
use crate::structures::Observation;
use crate::total_order::structures;
use std::cmp::Ordering;
use std::iter::{once, repeat_n};

pub fn algorithm<
    D: Direction,
    F: CauchyMeanValueFunctional,
    I: Into<Observation<f64, F::Response, F::Censoring>>,
>(
    data: impl Iterator<Item = I>,
    functional: &F,
) -> Vec<f64> {
    let allocated = structures::allocate_and_sort(data);
    algorithm_pre_sorted::<D, _, _, _>(allocated.into_iter(), functional)
}

/// Vanilla pool-adjacent-violators over covariate groups, generic in the functional. The
/// functional owns all numerical state: each block carries `F::Block`, and any algorithm-wide
/// scratch lives in `F::Shared`. The algorithm itself only manages the stack and the
/// ordering check.
pub fn algorithm_pre_sorted<
    D: Direction,
    F: CauchyMeanValueFunctional,
    Cov: PartialEq,
    I: Into<Observation<Cov, F::Response, F::Censoring>>,
>(
    iter: impl Iterator<Item = I>,
    functional: &F,
) -> Vec<f64> {
    // Pass 1: collect observations, compress by equal-covariate groups (covariates discarded).
    let mut data: Vec<Observation<(), F::Response, F::Censoring>> = Vec::new();
    let mut group_ends: Vec<usize> = Vec::new();
    let mut iter = iter.map(I::into);
    let mut next = iter.next();
    assert!(next.is_some(), "at least one item needed");

    while let Some(first) = next {
        let group_cov = first.x;
        data.push(Observation {
            x: (),
            y: first.y,
            observed: first.observed,
            weight: first.weight,
        });
        next = loop {
            match iter.next() {
                Some(o) if o.x == group_cov => {
                    data.push(Observation {
                        x: (),
                        y: o.y,
                        observed: o.observed,
                        weight: o.weight,
                    });
                }
                other => break other,
            }
        };
        group_ends.push(data.len());
    }

    let n_groups = group_ends.len();
    let get_data = |g: usize| {
        let start = if g > 0 { group_ends[g - 1] } else { 0 };
        let end = group_ends[g];
        data[start..end].iter().copied()
    };

    // Pass 2: vanilla PAVA. Stack entries are `(block_start, block)`, where the block's range
    // is `[block_start, next_block_start)` (or `[block_start, n_groups)` for the top of stack).
    let mut shared = functional.init_shared(n_groups);
    let mut stack: Vec<(usize, F::Block)> = Vec::with_capacity(n_groups);

    for g in 0..n_groups {
        stack.push((g, functional.singleton_block(g, &get_data, &mut shared)));
        // The newly pushed block always covers `[g, g + 1)`. Each merge below extends its
        // *start* leftward to absorb the popped neighbour, but the end stays at `g + 1`.
        let last_end = g + 1;

        while stack.len() >= 2 {
            let last_idx = stack.len() - 1;
            let before_last_idx = last_idx - 1;
            let bl_value = functional.block_value(&stack[before_last_idx].1);
            let last_value = functional.block_value(&stack[last_idx].1);
            if bl_value.partial_cmp(&last_value).unwrap() != D::FORBIDDEN_ORDERING {
                break;
            }
            let (last_start, last_block) = stack.pop().unwrap();
            let bl_idx = stack.len() - 1;
            let bl_start = stack[bl_idx].0;
            functional.merge_blocks(
                &mut stack[bl_idx].1,
                bl_start..last_start,
                last_block,
                last_start..last_end,
                &get_data,
                &mut shared,
            );
        }
    }

    // Output: one value per group, broadcast over each surviving block.
    let mut output = Vec::with_capacity(n_groups);
    for i in 0..stack.len() {
        let start = stack[i].0;
        let end = if i + 1 < stack.len() {
            stack[i + 1].0
        } else {
            n_groups
        };
        let value = functional.block_value(&stack[i].1);
        output.extend(repeat_n(value, end - start));
    }
    output
}

pub fn algorithm_definition<
    D: Direction,
    F: Functional,
    I: Into<Observation<f64, F::Response, F::Censoring>>,
>(
    data: impl Iterator<Item = I>,
    functional: &F,
) -> Vec<f64> {
    let allocated = structures::allocate_and_sort(data);
    algorithm_pre_sorted_definition::<D, _, _>(allocated, functional)
}

/// Uses the max-min definition
pub fn algorithm_pre_sorted_definition<
    D: Direction,
    F: Functional,
    Cov: Clone + Copy + PartialEq,
>(
    data: Vec<Observation<Cov, F::Response, F::Censoring>>,
    functional: &F,
) -> Vec<f64> {
    assert!(!data.is_empty());

    let start_indices: Vec<_> = data
        .iter()
        .enumerate()
        .filter(|&(i, item)| i == 0 || item.x != data[i - 1].x)
        .map(|(i, _)| i)
        .collect();

    let mut triangle = (0..start_indices.len())
        .map(|r| {
            (r..start_indices.len())
                .map(|s| {
                    let final_index = if s == start_indices.len() - 1 {
                        data.len()
                    } else {
                        start_indices[s + 1]
                    };
                    functional.evaluate(data[start_indices[r]..final_index].iter().copied())
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    for diagonal in 1..start_indices.len() {
        for diagonal_item in 0..start_indices.len() - diagonal {
            let r = diagonal_item;
            let s = diagonal + diagonal_item;

            for boundary in r..s {
                let left = triangle[r][boundary - r];
                let right = triangle[boundary + 1][s - (boundary + 1)];

                if left.is_nan() || right.is_nan() {
                    continue;
                }

                let (lowest, highest) = (left.min(right), left.max(right));

                triangle[r][s - r] = triangle[r][s - r].min(highest).max(lowest);
            }
        }
    }

    type MinOrMax = fn(f64, f64) -> f64;
    let (outer_direction, inner_direction): (MinOrMax, MinOrMax) = match D::FORBIDDEN_ORDERING {
        Ordering::Less => (f64::min, f64::max),
        Ordering::Equal => panic!(),
        Ordering::Greater => (f64::max, f64::min),
    };

    let extremes = (0..start_indices.len())
        .map(|i| {
            (0..=i)
                .map(|r| {
                    (i..start_indices.len())
                        .map(|s| triangle[r][s - r])
                        .reduce(inner_direction)
                        .unwrap()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    extremes
        .into_iter()
        .map(|r_values| r_values.into_iter().reduce(outer_direction).unwrap())
        .collect()
}

// Specialized implementations to analyze behavior

/// Compute the global variance. Used only for assessing convergence on an iid sample.
#[must_use]
pub fn algorithm_pre_sorted_global_variance(responses: &[f64]) -> f64 {
    let cumulative_sum: Vec<_> = responses
        .iter()
        .scan((0.0, 0.0), |(total, total_sq), value| {
            *total += value;
            *total_sq += value * value;
            Some((*total, *total_sq))
        })
        .collect();

    let mut triangle = (0..responses.len() - 1)
        .map(|r| {
            (r + 1..responses.len())
                .map(|s| {
                    let (total, total_sq) = if r > 0 {
                        (
                            cumulative_sum[s].0 - cumulative_sum[r - 1].0,
                            cumulative_sum[s].1 - cumulative_sum[r - 1].1,
                        )
                    } else {
                        (cumulative_sum[s].0, cumulative_sum[s].1)
                    };
                    let count = s - r + 1;

                    (total_sq - total * total / count as f64) / (count - 1) as f64
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    for diagonal in 4..=responses.len() {
        // diagonal is the interval length, 1..responses.len()

        for diagonal_item in 0..=responses.len() - diagonal {
            // diagonal item is the start index of the interval being modified

            // start index, inclusive
            let r = diagonal_item;
            // end index, inclusive
            let s = r + diagonal - 1;

            let median = median(
                ((r + 1)..(s - 1))
                    .flat_map(|k| [triangle[r][k - r - 1], triangle[k + 1][s - (k + 1) - 1]])
                    .chain(once(triangle[r][s - r - 1])),
            );
            triangle[r][s - r - 1] = median;
        }
    }

    *triangle[0].last().unwrap()
}

#[derive(Debug, PartialEq)]
pub struct VarianceSelection {
    pub value: f64,
    /// Inclusive start index in `responses` of the *leaf interval* whose (sample) variance
    /// ultimately determines `value` after the DP/median propagation.
    pub start: usize,
    /// Inclusive end index in `responses` of the *leaf interval* (see `start`).
    pub end: usize,
    /// Number of leaf blocks in the finest (binary-split) partition along the chosen trace.
    pub blocks: usize,
    /// Largest leaf-block length (in number of points) in that partition.
    pub largest: usize,
    /// Smallest leaf-block length (in number of points) in that partition.
    pub smallest: usize,
}

/// Compact per-cell trace metadata to reduce memory vs storing `usize` everywhere.
/// (You can switch these to `usize` if you really need `n > u32::MAX`.)
#[derive(Clone, Copy, Debug)]
struct TraceCell {
    start: u32,    // leaf interval start
    end: u32,      // leaf interval end
    blocks: u32,   // number of leaf blocks in traced partition
    largest: u32,  // max leaf-block length
    smallest: u32, // min leaf-block length
}

#[inline]
fn interval_len_u32(r: usize, s: usize) -> u32 {
    (s - r + 1) as u32
}

#[inline]
fn sample_variance_from_prefix(prefix: &[(f64, f64)], r: usize, s: usize) -> f64 {
    let (total, total_sq) = if r > 0 {
        (prefix[s].0 - prefix[r - 1].0, prefix[s].1 - prefix[r - 1].1)
    } else {
        (prefix[s].0, prefix[s].1)
    };

    let count = (s - r + 1) as f64;
    // unbiased sample variance
    (total_sq - total * total / count) / (count - 1.0)
}

#[must_use]
pub fn algorithm_pre_sorted_global_variance_with_trace(responses: &[f64]) -> VarianceSelection {
    let n = responses.len();
    assert!(n >= 2, "need at least 2 responses");
    assert!(
        n <= u32::MAX as usize,
        "n too large for u32 trace; switch TraceCell fields to usize if needed"
    );

    // Prefix sums: (sum, sumsq)
    let mut prefix = Vec::with_capacity(n);
    let mut total = 0.0;
    let mut total_sq = 0.0;
    for &x in responses {
        total += x;
        total_sq += x * x;
        prefix.push((total, total_sq));
    }

    // Triangle storage
    // row r stores intervals [r, s] for s=r+1..n-1 at index (s-r-1)
    let mut tri_val: Vec<Vec<f64>> = (0..n - 1)
        .map(|r| {
            let mut row = Vec::with_capacity(n - r - 1);
            for s in (r + 1)..n {
                row.push(sample_variance_from_prefix(&prefix, r, s));
            }
            row
        })
        .collect();

    let mut tri_tr: Vec<Vec<TraceCell>> = (0..n - 1)
        .map(|r| {
            let mut row = Vec::with_capacity(n - r - 1);
            for s in (r + 1)..n {
                let len = interval_len_u32(r, s);
                row.push(TraceCell {
                    start: r as u32,
                    end: s as u32,
                    blocks: 1,
                    largest: len,
                    smallest: len,
                });
            }
            row
        })
        .collect();

    for diagonal in 4..=n {
        // diagonal = interval length
        for r in 0..=n - diagonal {
            let s = r + diagonal - 1;
            let idx = s - r - 1;

            // Copy out active before overwrite (avoids borrow issues and keeps logic clear).
            let active_val = tri_val[r][idx];
            let active_tr = tri_tr[r][idx];

            let med = median(
                ((r + 1)..(s - 1))
                    .flat_map(|k| [tri_val[r][k - r - 1], tri_val[k + 1][s - (k + 1) - 1]])
                    .chain(once(active_val)),
            );

            // Choose a deterministic trace corresponding to one element that equals `med`.
            // (Median returns one of the candidate values; ties are possible, trace is then arbitrary.)
            let med_bits = med.to_bits();
            let matches = |x: f64| x.to_bits() == med_bits || x == med;

            let chosen_tr: TraceCell = if matches(active_val) {
                // Median selected S_[r:s]^2 (the interval's own sample variance).
                active_tr
            } else {
                // Otherwise, it must match one of the left-anchored or right-anchored candidates.
                let mut found: Option<TraceCell> = None;

                // Prefer left candidates first (deterministic tie-break).
                for k in (r + 1)..(s - 1) {
                    let left_val = tri_val[r][k - r - 1];
                    if matches(left_val) {
                        let lt = tri_tr[r][k - r - 1];
                        let right_len = (s - k) as u32; // length of [k+1, s]
                        found = Some(TraceCell {
                            start: lt.start,
                            end: lt.end,
                            blocks: lt.blocks + 1,
                            largest: lt.largest.max(right_len),
                            smallest: lt.smallest.min(right_len),
                        });
                        break;
                    }
                }

                // Then right candidates.
                if found.is_none() {
                    for k in (r + 1)..(s - 1) {
                        let rr = k + 1;
                        let right_val = tri_val[rr][s - rr - 1];
                        if matches(right_val) {
                            let rt = tri_tr[rr][s - rr - 1];
                            let left_len = (rr - r) as u32; // length of [r, k]
                            found = Some(TraceCell {
                                start: rt.start,
                                end: rt.end,
                                blocks: rt.blocks + 1,
                                largest: rt.largest.max(left_len),
                                smallest: rt.smallest.min(left_len),
                            });
                            break;
                        }
                    }
                }

                found.expect("median did not match any candidate (unexpected)")
            };

            // Store updated cell.
            tri_val[r][idx] = med;
            tri_tr[r][idx] = chosen_tr;
        }
    }

    let value = *tri_val[0].last().unwrap();
    let tr = *tri_tr[0].last().unwrap();

    VarianceSelection {
        value,
        start: tr.start as usize,
        end: tr.end as usize,
        blocks: tr.blocks as usize,
        largest: tr.largest as usize,
        smallest: tr.smallest as usize,
    }
}

/// Compute the global variance. Used only for assessing convergence on an iid sample.
#[must_use]
pub fn algorithm_pre_sorted_global_kaplan_meier(
    threshold: f64,
    responses: &[f64],
    censoring: &[bool],
) -> f64 {
    assert_eq!(responses.len(), censoring.len());
    let n = responses.len();
    let data: Vec<_> = responses
        .iter()
        .zip(censoring.iter())
        .map(|(&r, &c)| Observation {
            x: (),
            y: r,
            observed: c,
            weight: 0.0,
        })
        .collect();

    let functional = |slice: &[Observation<_, f64, bool>]| {
        let mut total_weight = slice.len();
        let mut subset: Vec<_> = slice.iter().filter(|o| o.y <= threshold).copied().collect();

        subset.sort_unstable_by(|l, r| {
            l.y.total_cmp(&r.y)
                .then(l.observed.cmp(&r.observed).reverse())
        });

        let mut survival = 1.0;
        for o in subset {
            if o.observed {
                survival *= 1.0 - 1.0 / total_weight as f64;
            }
            total_weight -= 1;
        }

        1.0 - survival
    };

    let mut triangle = (0..n)
        .map(|r| (r..n).map(|s| functional(&data[r..=s])).collect::<Vec<_>>())
        .collect::<Vec<_>>();

    for diagonal in 1..n {
        for diagonal_item in 0..n - diagonal {
            let r = diagonal_item;
            let s = diagonal + diagonal_item;

            for boundary in r..s {
                let left = triangle[r][boundary - r];
                let right = triangle[boundary + 1][s - (boundary + 1)];

                if left.is_nan() || right.is_nan() {
                    continue;
                }

                let (lowest, highest) = (left.min(right), left.max(right));

                triangle[r][s - r] = triangle[r][s - r].min(highest).max(lowest);
            }
        }
    }

    *triangle[0].last().unwrap()
}

#[cfg(test)]
mod test {
    use crate::functionals::{ClippingWrapper, ExpectedShortfall, KaplanMeier, Variance};
    use crate::structures::{Decreasing, Increasing};
    use crate::test::is_relative_eq_vec;
    use crate::total_order::functionals::{
        VarianceSelection, algorithm_definition, algorithm_pre_sorted,
        algorithm_pre_sorted_global_variance, algorithm_pre_sorted_global_variance_with_trace,
    };
    use approx::assert_relative_eq;
    use itertools::izip;
    use std::panic;

    #[test]
    fn test_variance_trivial() {
        let functional = Variance::new();
        assert_eq!(
            algorithm_definition::<Increasing, _, _>(
                [(1.0, 1.0), (2.0, 3.0)].into_iter(),
                &functional,
            ),
            vec![2.0, 2.0],
        )
    }

    #[test]
    fn test_variance_two_groups_ordering_satisfied() {
        let functional = Variance::new();
        let low = (0..=10).map(|i| (1.0, 100.0 + i as f64 / 10.0));
        let high = (0..=10).map(|i| (2.0, i as f64));

        let result = algorithm_definition::<Increasing, _, _>(low.chain(high), &functional);
        let expected = vec![0.11, 11.0];

        compare_arrays(&result, &expected);
    }

    #[test]
    fn test_variance_two_groups_ordering_violated() {
        let functional = Variance::new();
        let basis = (0..=10).map(|i| (1.0, i as f64));
        let shifted = (0..=10).map(|i| (2.0, 100.0 + i as f64 / 10.0));

        let result = algorithm_definition::<Increasing, _, _>(basis.chain(shifted), &functional);
        let expected = vec![11.0, 11.0];

        compare_arrays(&result, &expected);
    }

    #[test]
    fn test_expected_shortfall_two_groups_ordering_satisfied() {
        let functional = ExpectedShortfall::new(0.2);
        let low = (0..10).map(|i| (1.0, i as f64));
        let high = (0..10).map(|i| (2.0, 100.0 + i as f64));

        let result = algorithm_definition::<Increasing, _, _>(low.chain(high), &functional);
        let expected = vec![0.5, 100.0 + 0.5];

        compare_arrays(&result, &expected);
    }

    #[test]
    fn test_expected_shortfall_two_groups_ordering_violated() {
        let functional = ExpectedShortfall::new(0.2);
        let low = (0..10).map(|i| (1.0, 100.0 + i as f64));
        let high = (0..10).map(|i| (2.0, i as f64));

        let result = algorithm_definition::<Increasing, _, _>(low.chain(high), &functional);
        let expected = vec![1.5, 1.5];

        compare_arrays(&result, &expected);
    }

    fn compare_arrays(result: &[f64], expected: &[f64]) {
        assert_eq!(result.len(), expected.len());
        for (r, e) in result.into_iter().zip(expected.into_iter()) {
            assert_relative_eq!(r, e, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_global_variance() {
        let responses = vec![-1.0, 1.0, -1.0, 1.0];

        let result = algorithm_pre_sorted_global_variance(&responses);

        assert_relative_eq!(result, 2.0);
    }

    #[test]
    fn test_global_variance_uncapped() {
        let responses = vec![0.0, 0.0, 100.0, 0.0, 0.0];

        let result = algorithm_pre_sorted_global_variance(&responses);

        assert_relative_eq!(result, 2000.0);
    }

    #[test]
    fn test_algorithm_kaplan_meier() {
        assert_eq!(
            algorithm_pre_sorted::<Decreasing, _, _, _>(
                izip!(
                    [1.0, 2.0, 3.0, 4.0],
                    [2.0, 1.0, 4.0, 3.0],
                    [false, true, true, true],
                ),
                &ClippingWrapper::new(KaplanMeier::new(3.5)),
            ),
            vec![0.5, 0.5, 0.5, 0.5],
        );
        assert_eq!(
            algorithm_pre_sorted::<Decreasing, _, _, _>(
                izip!(
                    [1.0, 2.0, 3.0, 4.0],
                    [2.0, 1.0, 4.0, 3.0],
                    [false, true, true, true],
                ),
                &ClippingWrapper::new(KaplanMeier::new(2.5)),
            ),
            vec![0.5, 0.5, 0.0, 0.0],
        );
    }

    #[test]
    fn test_algorithm_kaplan_meier_2() {
        let result = algorithm_pre_sorted::<Decreasing, _, _, _>(
            izip!(
                [
                    1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
                    15.0, 16.0, 17.0, 18.0
                ],
                [
                    16.0_f64, 9.0, 17.0, 7.0, 10.0, 8.0, 2.0, 13.0, 15.0, 1.0, 14.0, 3.0, 11.0,
                    12.0, 4.0, 18.0, 6.0, 5.0
                ],
                [
                    false, false, false, true, false, true, true, false, false, true, true, true,
                    false, true, false, true, true, true
                ],
            ),
            &ClippingWrapper::new(KaplanMeier::new(13.0)),
        );
        let expected = [24.0 / 49.0; 18];
        assert!(
            is_relative_eq_vec(&result, &expected),
            "Result:   {:?}\nExpected: {:?}\n",
            result,
            expected,
        );
    }

    #[test]
    fn test_variance_with_trace() {
        let result = algorithm_pre_sorted_global_variance_with_trace(&[3.0, 3.0, 4.0, 4.00000001]);
        assert_eq!(
            result,
            VarianceSelection {
                value: 0.0,
                start: 0,
                end: 1,
                blocks: 2,
                largest: 2,
                smallest: 2,
            },
        );
        let result = algorithm_pre_sorted_global_variance_with_trace(&[1., -1., 1., -1. - 1e-2]);
        assert_eq!(
            result,
            VarianceSelection {
                value: 2.,
                start: 0,
                end: 1,
                blocks: 2,
                largest: 2,
                smallest: 2,
            },
        );
        let result = algorithm_pre_sorted_global_variance_with_trace(&[1., -1., 1., -1. + 1e-2]);
        assert_eq!(
            result,
            VarianceSelection {
                value: 1.98005,
                start: 2,
                end: 3,
                blocks: 2,
                largest: 2,
                smallest: 2,
            },
        );
        let result = algorithm_pre_sorted_global_variance_with_trace(&[1., -1., 0., 0.]);
        assert_eq!(
            result,
            VarianceSelection {
                value: 2. / 3.,
                start: 0,
                end: 3,
                blocks: 1,
                largest: 4,
                smallest: 4,
            },
        );
    }

    #[test]
    fn test_compare() {
        let data = [
            -0.83817967,
            0.20491239,
            -0.98949905,
            0.74742851,
            -0.26993208,
        ];
        let result = algorithm_pre_sorted_global_variance(&data);
        let with_trace = algorithm_pre_sorted_global_variance_with_trace(&data);
        assert!(
            (result - with_trace.value).abs() < 1e-3,
            "{} {}",
            result,
            with_trace.value
        );
    }
}
