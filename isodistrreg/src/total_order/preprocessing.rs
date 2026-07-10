use crate::Float;
use crate::error::Error;
use crate::preprocessing::validate;
use crate::structures::Observation;
use crate::total_order::structures::{AlgorithmContext, CensoredContext, CovariateStatistic};
use itertools::{Itertools, izip};

pub fn preprocess_uncensored<X: Float, Y: Float, W: Float>(
    x: &[X],
    y: &[Y],
    weight: &[W],
) -> AlgorithmContext<X, Y, ()> {
    let mut context = preprocess(x, y, |_| (), weight);
    context.unique_responses = context
        .observations
        .first()
        .map(|o| o.y)
        .into_iter()
        .chain(
            context
                .observations
                .array_windows()
                .filter_map(|[left, right]| (right.y != left.y).then_some(right.y)),
        )
        .collect();
    context
}

/// Finalize the in-progress (last) observation: store its narrowed weight and fold it into
/// the current (last) covariate statistic.
///
/// A positive `W`-precision weight sum can narrow to exactly 0.0f32 (a `W = f64` weight
/// below half the smallest f32 subnormal). Such an observation carries no representable
/// mass in the f32 algorithm domain and is dropped here exactly like a zero-weight
/// observation (see the `fit()` weight contract), keeping the "every weight in the context
/// is positive" invariant intact.
fn commit_observation<Y, S>(
    observations: &mut Vec<Observation<usize, Y, S, f32>>,
    covariate_statistics: &mut [CovariateStatistic],
    finalized: f32,
) {
    if finalized > 0.0 {
        observations.last_mut().unwrap().weight = finalized;
        let statistic = covariate_statistics.last_mut().unwrap();
        statistic.weight += finalized;
        statistic.cumulative_weight += finalized;
    } else {
        observations.pop();
    }
}

/// Close the covariate whose observations have all been committed: if every one of its
/// observations was dropped (narrowed to 0.0f32), the covariate carries no mass and must not
/// enter the grid — exactly like a covariate whose observations all have zero weight.
fn close_covariate<X>(
    covariate_statistics: &mut Vec<CovariateStatistic>,
    unique_covariates: &mut Vec<X>,
) {
    if covariate_statistics.last().unwrap().weight == 0.0 {
        covariate_statistics.pop();
        unique_covariates.pop();
    }
}

/// Build the algorithm context: observations sorted and deduplicated, covariate grid and
/// per-covariate weight statistics.
///
/// Zero-weight observations are dropped here, before anything is aggregated: they carry no
/// statistical information (every estimator the kernels compute is invariant to adding
/// them), their response values create no thresholds, and a covariate whose observations
/// all have zero weight does not enter the grid. Downstream kernels divide by observation
/// and covariate weights and rely on every weight in the context being positive. With all
/// observations dropped, the returned context is empty (the empty fit).
///
/// Weights are narrowed to f32 here (see the `fit()` weight contract); positive weights
/// that narrow to 0.0f32 are dropped like zero weights ([`commit_observation`]), with the
/// covariate grid and threshold grid matching the post-drop subsample exactly.
pub fn preprocess<X: Float, Y: Float, S: Ord, W: Float, F: Fn(usize) -> S>(
    x: &[X],
    y: &[Y],
    observed: F,
    weight: &[W],
) -> AlgorithmContext<X, Y, S> {
    let n = x.len();

    // Sort packed values rather than argsort indices, for the same reason as in
    // `preprocess_censored`: the comparator reads its keys from the element itself
    // instead of gathering from separate arrays per comparison, and the aggregation
    // loop below scans one contiguous buffer instead of gathering through the sorted
    // order. Same comparator, still unstable — order among fully-tied keys remains
    // arbitrary. `observed` is evaluated once per row at pack time (it is a pure
    // function of the original index, which filtering does not disturb).
    //
    // Zero-weight observations are dropped while packing, before they are even
    // compared. `validate` has rejected negative and non-finite weights, so `> 0` is
    // exactly "nonzero" here.
    struct Item<X, Y, S, W> {
        x: X,
        y: Y,
        observed: S,
        weight: W,
    }
    // Canonicalize signed zeros while packing (`-0.0 + 0.0 == +0.0` in IEEE
    // round-to-nearest; every other value is unchanged): the sort keys use `total_cmp`,
    // which separates -0.0 from +0.0, while the aggregation below merges them with `==`.
    // Mixed zeros would otherwise assemble one covariate or threshold out of two
    // sort-distinct groups — leaving duplicate observations non-adjacent (escaping
    // deduplication) and a threshold's events not ascending in covariate, both contracts
    // the PAVA batching depends on.
    let mut items: Vec<Item<X, Y, S, W>> = izip!(x, y, weight)
        .enumerate()
        .filter(|&(_, (_, _, &weight))| weight > W::zero())
        .map(|(index, (&x, &y, &weight))| Item {
            x: x + X::zero(),
            y: y + Y::zero(),
            observed: observed(index),
            weight,
        })
        .collect();
    if items.is_empty() {
        return AlgorithmContext {
            observations: Vec::with_capacity(0),
            covariate_statistics: Vec::with_capacity(0),
            unique_responses: Vec::with_capacity(0),
            unique_covariates: Vec::with_capacity(0),
        };
    }
    items.sort_unstable_by(|a, b| {
        a.x.total_cmp(&b.x)
            .then(a.y.total_cmp(&b.y))
            .then(a.observed.cmp(&b.observed).reverse())
    });

    // While copying the data in the sorted order and aggregating identical observations by weight,
    // we track the unique covariates we encounter.
    let mut observations: Vec<Observation<usize, Y, S, f32>> = Vec::with_capacity(n);
    let mut unique_covariates = Vec::with_capacity(n);
    let mut covariate_statistics: Vec<CovariateStatistic> = Vec::with_capacity(n);

    // Accumulator for the in-progress observation's weight in the caller's W precision.
    // Invariant: `last_w_accum` holds the running W-precision sum that will be committed
    // into `observations.last().weight` (narrowed to f32) when the in-progress observation
    // is finalized (i.e., when a new observation is pushed or at the end of the loop).
    // While the in-progress observation is current, its `weight: f32` field holds an
    // uninitialized placeholder (0.0) and must not be read.
    let mut last_w_accum: W;

    let mut items_sorted = items.into_iter();
    {
        let item = items_sorted.next().unwrap();
        observations.push(Observation {
            x: 0,
            y: item.y,
            observed: item.observed,
            weight: 0.0, // placeholder; finalized from `last_w_accum`
        });
        last_w_accum = item.weight;
        covariate_statistics.push(CovariateStatistic {
            weight: 0.0,
            cumulative_weight: 0.0, // Accumulated as each observation's weight is committed
        });
        unique_covariates.push(item.x);
    }

    for item in items_sorted {
        let covariate_equal = item.x == *unique_covariates.last().unwrap();
        let (response_equal, censoring_equal) = {
            let last_observation = observations.last().unwrap();
            (
                item.y == last_observation.y,
                item.observed == last_observation.observed,
            )
        };

        if covariate_equal && response_equal && censoring_equal {
            // At the same observation -> just accumulate observation weight in W precision
            last_w_accum = last_w_accum + item.weight;
        } else if covariate_equal {
            // New observation but same covariate -> finalize the previous observation's
            // weight (narrow to f32, dropping it if no mass survives the narrowing),
            // accumulate covariate weight from it, and start a new observation with a
            // fresh W accumulator.
            commit_observation(
                &mut observations,
                &mut covariate_statistics,
                last_w_accum.to_f32().unwrap(),
            );
            observations.push(Observation {
                x: unique_covariates.len() - 1, // we stay with the same covariate
                y: item.y,
                observed: item.observed,
                weight: 0.0, // placeholder; finalized from `last_w_accum`
            });
            last_w_accum = item.weight;
        } else {
            // New observation and new covariate -> finalize the previous observation's
            // weight, accumulate covariate weight from it, close the finished covariate
            // (dropping it if all its observations were dropped), push a new covariate
            // statistic chaining the cumulative weight, push a new observation, and seed
            // a fresh W accumulator.
            commit_observation(
                &mut observations,
                &mut covariate_statistics,
                last_w_accum.to_f32().unwrap(),
            );
            close_covariate(&mut covariate_statistics, &mut unique_covariates);
            let new_statistic = CovariateStatistic {
                weight: 0.0,
                cumulative_weight: covariate_statistics
                    .last()
                    .map_or(0.0, |statistic| statistic.cumulative_weight),
            };
            covariate_statistics.push(new_statistic);
            observations.push(Observation {
                x: unique_covariates.len(),
                y: item.y,
                observed: item.observed,
                weight: 0.0, // placeholder; finalized from `last_w_accum`
            });
            last_w_accum = item.weight;
            unique_covariates.push(item.x);
        }
    }
    // Finalize the last in-progress observation and close the last covariate.
    commit_observation(
        &mut observations,
        &mut covariate_statistics,
        last_w_accum.to_f32().unwrap(),
    );
    close_covariate(&mut covariate_statistics, &mut unique_covariates);
    // Could shrink to fit `observations`, `unique_covariates` and `covariate_statistics` here if
    // duplicates are common, and we care about memory usage

    // Stable to keep each response threshold sorted by covariate - this is probably faster in PAVA
    observations.sort_by(|l, r| {
        l.y.total_cmp(&r.y)
            .then(l.observed.cmp(&r.observed).reverse())
    });

    AlgorithmContext {
        observations,
        covariate_statistics,
        unique_responses: Vec::with_capacity(0),
        unique_covariates,
    }
}

/// Preprocess censored data such that the censored total-order algorithms can be run.
/// Shared by the stochastic-dominance and hazard-rate pipelines.
///
/// Thresholds contain all unique response values that have at least one uncensored
/// observation. Censored observations sorting before the first event carry no information
/// at any threshold (they leave the risk set before every event) and are discarded —
/// covariates whose mass is entirely discarded this way drop out of the grid. Each
/// remaining censored observation is snapped to the latest threshold at or below its
/// response, which leaves every interval Kaplan-Meier estimator unchanged (with the
/// events-before-censorings tie convention, censored mass at a threshold stays at risk
/// through that threshold's events).
///
/// Censored observations at the same covariate with no event between their responses —
/// equivalently, with the same snapped threshold — are merged into one observation per
/// (threshold, covariate), regardless of how other covariates' observations interleave
/// with them in the response order. The merge is exact for every estimator downstream: a
/// censored observation never carries a value, it only leaves the at-risk set, and the
/// at-risk weight is read only at event folds — so all of a threshold gap's censored
/// observations act between the same two events, where equal covariates make them
/// indistinguishable from their combined weight. (The hazard-rate pipeline folds censored
/// observations multiplicatively per covariate, where the same-covariate factors
/// telescope: (1−w₁/R)·(1−w₂/(R−w₁)) = 1−(w₁+w₂)/R.) A censored pair with an event
/// between them snaps to different thresholds and is never merged. Together with the
/// per-(threshold, covariate) event dedup this bounds the observation stream by 2·t·m,
/// which is what keeps the fits' per-cell walk work bounded by the pool sweeps' O(t·m³)
/// instead of scaling with the raw censored count.
pub fn preprocess_censored<X: Float, Y: Float, W: Float>(
    x: &[X],
    y: &[Y],
    observed: &[bool],
    weights: &[W],
) -> Result<CensoredContext<X, Y>, Error> {
    let n = validate(x.chunks_exact(1), y, Some(observed), Some(weights))?;

    let (observations_response_sorted, thresholds) = {
        // Sort packed values rather than argsort indices: the comparator then reads its
        // keys from the element itself instead of three gathers per comparison, which is
        // what made the index sort dominate preprocessing for large inputs. Same
        // comparator either way, both unstable — order among fully-tied keys remains
        // arbitrary.
        //
        // Zero-weight observations are dropped before anything is aggregated (or even
        // sorted): they carry no statistical information, create no thresholds, and the
        // kernels rely on every weight in the context being positive. `validate` has
        // rejected negative and non-finite weights, so `> 0` is exactly "nonzero".
        struct Item<X, Y, W> {
            y: Y,
            x: X,
            weight: W,
            observed: bool,
        }
        // Canonicalize signed zeros while packing, exactly as in `preprocess`: the
        // `total_cmp` sort keys must agree with the `==` aggregation below on ±0.0, or a
        // merged threshold's events end up not ascending in covariate and duplicate
        // covariates escape deduplication — contracts the batched PAVA updates rely on.
        let mut items: Vec<Item<X, Y, W>> = izip!(x, y, observed, weights)
            .filter(|&(_, _, _, &weight)| weight > W::zero())
            .map(|(&x, &y, &observed, &weight)| Item {
                y: y + Y::zero(),
                x: x + X::zero(),
                weight,
                observed,
            })
            .collect();
        items.sort_unstable_by(|a, b| {
            a.y.total_cmp(&b.y)
                .then(a.observed.cmp(&b.observed).reverse())
                .then(a.x.total_cmp(&b.x))
        });

        // Discard censored observations not greater than or equal to any uncensored observation
        let first_uncensored = items
            .iter()
            .find_position(|item| item.observed)
            .map(|(index, _)| index);
        let Some(first_uncensored_index) = first_uncensored else {
            return Ok(CensoredContext {
                observations: Vec::with_capacity(0),
                covariate_statistics: Vec::with_capacity(0),
                unique_covariates: Vec::with_capacity(0),
                thresholds: Vec::with_capacity(0),
            });
        };
        let capacity_upper_bound = n - first_uncensored_index;
        let mut obs = Vec::with_capacity(capacity_upper_bound);
        let mut thresholds = Vec::with_capacity(capacity_upper_bound);
        debug_assert!(items[first_uncensored_index].observed);

        // From the first event on, the stream alternates maximal runs of two kinds: one
        // threshold's events (equal response, covariates ascending, duplicates adjacent)
        // and one threshold gap's censored observations (everything until the next event,
        // sorted by raw response then covariate, so equal covariates need not be
        // adjacent). Both kinds aggregate per covariate, so one run = one batch: gather
        // (covariate, weight) pairs with weights in W precision, sort the censored batch
        // by covariate to make its interleaved duplicates adjacent, merge equal
        // covariates (narrowing each total to f32 once), and emit under the current
        // threshold index — which an event run has just pushed, and a censored run
        // inherits as its snapping target.
        let mut scratch: Vec<(X, W)> = Vec::new();
        for run in items[first_uncensored_index..]
            .chunk_by(|a, b| a.observed == b.observed && (!a.observed || a.y == b.y))
        {
            let observed = run[0].observed;
            if observed {
                thresholds.push(run[0].y);
            }
            for item in run {
                // Adjacent duplicates (all event duplicates, and censored duplicates
                // tied on the raw response) accumulate immediately.
                match scratch.last_mut() {
                    Some((x, weight)) if *x == item.x => *weight = *weight + item.weight,
                    _ => scratch.push((item.x, item.weight)),
                }
            }
            if !observed {
                scratch.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
            }
            let y = thresholds.len() - 1;
            let mut batch = scratch.drain(..);
            let (mut covariate, mut w_accum) = batch.next().unwrap();
            for (x, weight) in batch {
                if x == covariate {
                    w_accum = w_accum + weight;
                } else {
                    obs.push(Observation {
                        x: covariate,
                        y,
                        observed,
                        weight: w_accum.to_f32().unwrap(),
                    });
                    (covariate, w_accum) = (x, weight);
                }
            }
            obs.push(Observation {
                x: covariate,
                y,
                observed,
                weight: w_accum.to_f32().unwrap(),
            });
        }

        (obs, thresholds)
    };

    // Order of the observations by covariate (unstable is fine; we keep the response
    // order intact). Sorting (key, position) pairs beats an argsort for the same
    // gather-avoidance reason as above. Positions fit u32 — the estimator triangles
    // downstream could never be allocated for sizes anywhere near that limit.
    assert!(observations_response_sorted.len() <= u32::MAX as usize);
    let mut covariate_order: Vec<(X, u32)> = observations_response_sorted
        .iter()
        .enumerate()
        .map(|(position, o)| (o.x, position as u32))
        .collect();
    covariate_order.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

    let mut observations: Vec<Observation<usize, usize, bool, f32>> = vec![
        Observation {
            x: usize::MAX,
            y: usize::MAX,
            observed: false,
            weight: f32::NAN,
        };
        observations_response_sorted.len()
    ];
    let capacity_upper_bound = observations_response_sorted.len();
    let mut unique_covariates = Vec::with_capacity(capacity_upper_bound);
    let mut covariate_statistics = Vec::with_capacity(capacity_upper_bound);

    // First item
    let data_index = covariate_order[0].1 as usize;
    let observation = &observations_response_sorted[data_index];
    observations[data_index] = Observation {
        x: 0,
        y: observation.y,
        observed: observation.observed,
        weight: observation.weight,
    };
    unique_covariates.push(observation.x);
    covariate_statistics.push(CovariateStatistic {
        weight: observation.weight,
        cumulative_weight: 0.0,
    });
    // Remaining items
    for &(_, position) in covariate_order[1..].iter() {
        let data_index = position as usize;
        let observation = &observations_response_sorted[data_index];
        if observation.x != *unique_covariates.last().unwrap() {
            unique_covariates.push(observation.x);
            let last_statistic = covariate_statistics.last_mut().unwrap();
            last_statistic.cumulative_weight += last_statistic.weight;
            covariate_statistics.push(CovariateStatistic {
                weight: 0.0,
                cumulative_weight: covariate_statistics.last().unwrap().cumulative_weight,
            });
        }
        observations[data_index] = Observation {
            x: unique_covariates.len() - 1,
            y: observation.y,
            observed: observation.observed,
            weight: observation.weight,
        };
        covariate_statistics.last_mut().unwrap().weight += observation.weight;
    }
    let last_statistic = covariate_statistics.last_mut().unwrap();
    last_statistic.cumulative_weight += last_statistic.weight;
    debug_assert!(observations.iter().all(|o| !o.weight.is_nan()));

    Ok(CensoredContext {
        observations,
        covariate_statistics,
        unique_covariates,
        thresholds,
    })
}

#[cfg(test)]
mod test_censored {
    use crate::structures::Observation;
    use crate::total_order::preprocessing::preprocess_censored as preprocess;
    use crate::total_order::structures::{CensoredContext, CovariateStatistic};

    #[test]
    fn test_trivial_single_observation() {
        let context = preprocess(&[5.0], &[6.0], &[true], &[2.0]).ok().unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![Observation {
                    x: 0,
                    y: 0,
                    observed: true,
                    weight: 2.0,
                },],
                covariate_statistics: vec![CovariateStatistic {
                    weight: 2.0,
                    cumulative_weight: 2.0,
                },],
                unique_covariates: vec![5.0],
                thresholds: vec![6.0],
            }
        );
    }

    #[test]
    fn test_trivial_single_covariate() {
        let context = preprocess(&[5.0, 5.0], &[6.5, 6.0], &[true, true], &[1.0, 2.0])
            .ok()
            .unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 2.0,
                    },
                    Observation {
                        x: 0,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![CovariateStatistic {
                    weight: 3.0,
                    cumulative_weight: 3.0,
                },],
                unique_covariates: vec![5.0],
                thresholds: vec![6.0, 6.5],
            }
        );
    }

    #[test]
    fn test_trivial_single_response() {
        let context = preprocess(&[5.0, 7.0], &[6.0, 6.0], &[true, true], &[1.0, 3.0])
            .ok()
            .unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 1,
                        y: 0,
                        observed: true,
                        weight: 3.0,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 1.0,
                    },
                    CovariateStatistic {
                        weight: 3.0,
                        cumulative_weight: 4.0,
                    },
                ],
                unique_covariates: vec![5.0, 7.0],
                thresholds: vec![6.0],
            }
        );
    }

    #[test]
    fn test_monotone() {
        let context = preprocess(
            &[2.0, 1.0, 3.0],
            &[2.0, 1.0, 3.0],
            &[false, true, true],
            &[1.0, 9.0, 1.0],
        )
        .ok()
        .unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 9.0,
                    },
                    Observation {
                        x: 1,
                        y: 0,
                        observed: false,
                        weight: 1.0,
                    },
                    Observation {
                        x: 2,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 9.0,
                        cumulative_weight: 9.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 10.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 11.0,
                    },
                ],
                unique_covariates: vec![1.0, 2.0, 3.0],
                thresholds: vec![1.0, 3.0],
            }
        );
    }

    #[test]
    fn test_duplicate_response() {
        let context =
            preprocess(&[1.0, 3.0, 2.0], &[3.0, 3.0, 5.0], &[true; 3], &[1.0; 3]).unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 2,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 1,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 1.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 2.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 3.0,
                    },
                ],
                unique_covariates: vec![1.0, 2.0, 3.0],
                thresholds: vec![3.0, 5.0],
            },
        );
    }

    #[test]
    fn test_not_monotone_censored_2() {
        let context = preprocess(&[1.0, 2.0], &[4.0, 3.0], &[true, false], &[1.0, 1.0]).unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![Observation {
                    x: 0,
                    y: 0,
                    observed: true,
                    weight: 1.0,
                },],
                covariate_statistics: vec![CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 1.0,
                },],
                unique_covariates: vec![1.0],
                thresholds: vec![4.0],
            }
        );
    }

    #[test]
    fn test_censored_duplicates_merge_across_interleaving() {
        // Censored 2.3 and 2.7 at covariate 5.0 with a censored 2.5 at covariate 9.0
        // sorting between them: no event lies in (2.3, 2.7], so all three snap to the
        // threshold of the event at 2.0 and the covariate-5.0 pair merges despite not
        // being adjacent in the response order. Dyadic weights keep every sum exact.
        let context = preprocess(
            &[5.0, 5.0, 9.0, 5.0],
            &[2.0, 2.3, 2.5, 2.7],
            &[true, false, false, false],
            &[1.0, 0.25, 1.5, 0.5],
        )
        .unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 0,
                        y: 0,
                        observed: false,
                        weight: 0.75,
                    },
                    Observation {
                        x: 1,
                        y: 0,
                        observed: false,
                        weight: 1.5,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 1.75,
                        cumulative_weight: 1.75,
                    },
                    CovariateStatistic {
                        weight: 1.5,
                        cumulative_weight: 3.25,
                    },
                ],
                unique_covariates: vec![5.0, 9.0],
                thresholds: vec![2.0],
            }
        );
    }

    #[test]
    fn test_censored_duplicates_split_by_event() {
        // Same covariate pair, but now an event at 2.5 lies in (2.3, 2.7]: the two
        // censored observations snap to different thresholds — the 2.3 one leaves the
        // risk set before the 2.5 event, the 2.7 one after — and must NOT merge.
        let context = preprocess(
            &[5.0, 5.0, 9.0, 5.0],
            &[2.0, 2.3, 2.5, 2.7],
            &[true, false, true, false],
            &[1.0, 0.25, 1.5, 0.5],
        )
        .unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 0,
                        y: 0,
                        observed: false,
                        weight: 0.25,
                    },
                    Observation {
                        x: 1,
                        y: 1,
                        observed: true,
                        weight: 1.5,
                    },
                    Observation {
                        x: 0,
                        y: 1,
                        observed: false,
                        weight: 0.5,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 1.75,
                        cumulative_weight: 1.75,
                    },
                    CovariateStatistic {
                        weight: 1.5,
                        cumulative_weight: 3.25,
                    },
                ],
                unique_covariates: vec![5.0, 9.0],
                thresholds: vec![2.0, 2.5],
            }
        );
    }

    #[test]
    fn test_interleaved_and_adjacent_censored_duplicates_agree() {
        // The merge key is (snapped threshold, covariate): moving a censored response
        // within its threshold gap (2.7 -> 2.3 here, no event between either way) must
        // not change the preprocessed context at all.
        let interleaved = preprocess(
            &[5.0, 5.0, 9.0, 5.0],
            &[2.0, 2.3, 2.5, 2.7],
            &[true, false, false, false],
            &[1.0, 0.25, 1.5, 0.5],
        )
        .unwrap();
        let adjacent = preprocess(
            &[5.0, 5.0, 9.0, 5.0],
            &[2.0, 2.3, 2.5, 2.3],
            &[true, false, false, false],
            &[1.0, 0.25, 1.5, 0.5],
        )
        .unwrap();
        assert_eq!(interleaved, adjacent);
    }

    #[test]
    fn test_not_monotone_3() {
        let context = preprocess(
            &[2.0, 1.0, 3.0],
            &[3.0, 1.0, 2.0],
            &[false, true, true],
            &[1.0, 9.0, 1.0],
        )
        .unwrap();
        assert_eq!(
            context,
            CensoredContext {
                observations: vec![
                    Observation {
                        x: 0,
                        y: 0,
                        observed: true,
                        weight: 9.0,
                    },
                    Observation {
                        x: 2,
                        y: 1,
                        observed: true,
                        weight: 1.0,
                    },
                    Observation {
                        x: 1,
                        y: 1,
                        observed: false,
                        weight: 1.0,
                    },
                ],
                covariate_statistics: vec![
                    CovariateStatistic {
                        weight: 9.0,
                        cumulative_weight: 9.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 10.0,
                    },
                    CovariateStatistic {
                        weight: 1.0,
                        cumulative_weight: 11.0,
                    },
                ],
                unique_covariates: vec![1.0, 2.0, 3.0],
                thresholds: vec![1.0, 2.0],
            }
        );
    }
}

#[cfg(test)]
mod test {
    use crate::structures::Observation;
    use crate::total_order::preprocessing::preprocess_uncensored;
    use crate::total_order::structures::{AlgorithmContext, CovariateStatistic};

    #[test]
    fn test_trivial_single_observation() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[5.0], &[6.0], &[2.0]);
        assert_eq!(
            observations,
            vec![Observation {
                x: 0,
                y: 6.0,
                observed: (),
                weight: 2.0,
            },],
        );
        assert_eq!(
            covariate_statistics,
            vec![CovariateStatistic {
                weight: 2.0,
                cumulative_weight: 2.0,
            },],
        );
        assert_eq!(unique_responses, vec![6.0]);
        assert_eq!(unique_covariates, vec![5.0]);
    }

    #[test]
    fn test_trivial_single_covariate() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[5.0, 5.0], &[6.5, 6.0], &[1.0, 2.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 6.0,
                    observed: (),
                    weight: 2.0,
                },
                Observation {
                    x: 0,
                    y: 6.5,
                    observed: (),
                    weight: 1.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![CovariateStatistic {
                weight: 3.0,
                cumulative_weight: 3.0,
            }],
        );
        assert_eq!(unique_responses, vec![6.0, 6.5]);
        assert_eq!(unique_covariates, vec![5.0]);
    }

    #[test]
    fn test_trivial_single_response() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[5.0, 7.0], &[6.0, 6.0], &[1.0, 3.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 6.0,
                    observed: (),
                    weight: 1.0,
                },
                Observation {
                    x: 1,
                    y: 6.0,
                    observed: (),
                    weight: 3.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 1.0,
                },
                CovariateStatistic {
                    weight: 3.0,
                    cumulative_weight: 4.0,
                },
            ],
        );
        assert_eq!(unique_responses, vec![6.0]);
        assert_eq!(unique_covariates, vec![5.0, 7.0]);
    }

    #[test]
    fn test_monotone() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[2.0, 1.0, 3.0], &[2.0, 1.0, 3.0], &[1.0, 9.0, 1.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 1.0,
                    observed: (),
                    weight: 9.0,
                },
                Observation {
                    x: 1,
                    y: 2.0,
                    observed: (),
                    weight: 1.0,
                },
                Observation {
                    x: 2,
                    y: 3.0,
                    observed: (),
                    weight: 1.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![
                CovariateStatistic {
                    weight: 9.0,
                    cumulative_weight: 9.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 10.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 11.0,
                },
            ],
        );
        assert_eq!(unique_responses, vec![1.0, 2.0, 3.0]);
        assert_eq!(unique_covariates, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_not_monotone_3() {
        let AlgorithmContext {
            observations,
            covariate_statistics,
            unique_responses,
            unique_covariates,
        } = preprocess_uncensored(&[2.0, 1.0, 3.0], &[3.0, 1.0, 2.0], &[1.0, 9.0, 1.0]);
        assert_eq!(
            observations,
            vec![
                Observation {
                    x: 0,
                    y: 1.0,
                    observed: (),
                    weight: 9.0,
                },
                Observation {
                    x: 2,
                    y: 2.0,
                    observed: (),
                    weight: 1.0,
                },
                Observation {
                    x: 1,
                    y: 3.0,
                    observed: (),
                    weight: 1.0,
                },
            ],
        );
        assert_eq!(
            covariate_statistics,
            vec![
                CovariateStatistic {
                    weight: 9.0,
                    cumulative_weight: 9.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 10.0,
                },
                CovariateStatistic {
                    weight: 1.0,
                    cumulative_weight: 11.0,
                },
            ],
        );
        assert_eq!(unique_responses, vec![1.0, 2.0, 3.0]);
        assert_eq!(unique_covariates, vec![1.0, 2.0, 3.0]);
    }
}
