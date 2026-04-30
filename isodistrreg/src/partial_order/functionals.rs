use crate::functionals::CauchyMeanValueFunctional;
use crate::partial_order::{BitSet, routines};
use crate::structures::{Direction, Observation};
use std::collections::HashMap;
use std::iter::once;

/// Main algorithm: computes for each i:
///   max_{U upper, i∈U} min_{L lower, i∈L} M_{L∩U}
///
/// Uses only the oracle `mean` and pooling-betweenness pruning for the inner minimization.
pub fn algorithm_pre_sorted<
    D: Direction,
    F: CauchyMeanValueFunctional,
    I: Into<Observation<(), F::Response, F::Censoring>>,
>(
    observations: impl ExactSizeIterator<Item = I>,
    edges: &[(usize, usize)],
    functional: &F,
) -> Vec<f64> {
    match (observations.len() - 1) / BitSet::<1>::capacity() {
        0..1 => algorithm_pre_sorted_sized::<D, _, 1, _>(observations, edges, functional),
        1..2 => algorithm_pre_sorted_sized::<D, _, 2, _>(observations, edges, functional),
        2..4 => algorithm_pre_sorted_sized::<D, _, 4, _>(observations, edges, functional),
        4..8 => algorithm_pre_sorted_sized::<D, _, 8, _>(observations, edges, functional),
        8..16 => algorithm_pre_sorted_sized::<D, _, 16, _>(observations, edges, functional),
        16..32 => algorithm_pre_sorted_sized::<D, _, 32, _>(observations, edges, functional),
        32..64 => algorithm_pre_sorted_sized::<D, _, 64, _>(observations, edges, functional),
        64..128 => algorithm_pre_sorted_sized::<D, _, 128, _>(observations, edges, functional),
        _ => unimplemented!(
            "data set is too large, largest supported is n = {}",
            BitSet::<128>::capacity(),
        ),
    }
}

/// Requires unique covariates that are sorted
fn algorithm_pre_sorted_sized<
    D: Direction,
    F: CauchyMeanValueFunctional,
    const B: usize,
    I: Into<Observation<(), F::Response, F::Censoring>>,
>(
    observations: impl ExactSizeIterator<Item = I>,
    edges: &[(usize, usize)],
    functional: &F,
) -> Vec<f64>
where
    F::Response: Default,
{
    let n = observations.len();
    assert!(n <= BitSet::<B>::capacity());

    let (successors, predecessors, topological_order, _) =
        routines::compute_transitive_closure::<B>(n, edges);
    let mut topological_observations = vec![
        Observation {
            x: (),
            y: F::Response::default(),
            observed: F::Censoring::default(),
            weight: f64::NAN,
        };
        n
    ];
    for (i, observation) in observations.into_iter().enumerate() {
        topological_observations[topological_order[i]] = observation.into();
    }
    debug_assert!(topological_observations.iter().all(|o| !o.weight.is_nan()));

    // Inclusive closures for convenience
    let successors_inclusive = successors
        .into_iter()
        .enumerate()
        .map(|(i, s)| s.with_singleton(i))
        .collect::<Vec<_>>();

    algorithm_pre_sorted_inner::<D, _, _, _, _>(
        n,
        &topological_order,
        &successors_inclusive,
        &predecessors,
        &|i| once(topological_observations[i]),
        functional,
    )
}

pub fn algorithm_pre_sorted_inner<
    D: Direction,
    const B: usize,
    F: CauchyMeanValueFunctional,
    G: Fn(usize) -> I + Clone,
    I: Iterator<Item = Observation<(), F::Response, F::Censoring>> + Clone,
>(
    n: usize,
    topological_order: &[usize],
    successors_inclusive: &[BitSet<B>],
    predecessors: &[BitSet<B>],
    get_data: &G,
    functional: &F,
) -> Vec<f64> {
    // Store the sequence of lower sets and their values
    let mut out = vec![f64::NAN; n];
    let mut a_prev = if D::IS_INCREASING {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    };

    // Cache may or may not be used - depends on the functional
    let mut cache = HashMap::new();

    let mut active = BitSet::<B>::fill(n);
    while !active.is_empty() {
        let (chosen, best_value) = select_extreme_lower_set::<D, _, _, _, _>(
            &active,
            successors_inclusive,
            predecessors,
            get_data,
            functional,
            &mut cache,
        );

        a_prev = if D::IS_INCREASING {
            f64::max(best_value, a_prev)
        } else {
            f64::min(best_value, a_prev)
        };
        for item in chosen.iter() {
            out[topological_order[item]] = a_prev;
        }

        active = active.difference(&chosen);
    }
    debug_assert!(out.iter().all(|v| !v.is_nan()));
    out
}

/// Output-sensitive enumeration of all lower sets (ideals) of (V, ⪯).
fn select_extreme_lower_set<
    D: Direction,
    const B: usize,
    F: CauchyMeanValueFunctional,
    G: Fn(usize) -> I + Clone,
    I: Iterator<Item = Observation<(), F::Response, F::Censoring>> + Clone,
>(
    active: &BitSet<B>,
    successors_inclusive: &[BitSet<B>],
    predecessors: &[BitSet<B>],
    get_data: G,
    functional: &F,
    cache: &mut HashMap<BitSet<B>, f64>,
) -> (BitSet<B>, f64) {
    let mut best_value = if D::IS_INCREASING {
        f64::INFINITY
    } else {
        f64::NEG_INFINITY
    };
    let mut contenders = Vec::new();

    dfs::<D, _, _, _, _>(
        BitSet::new(),
        &BitSet::new(),
        &mut best_value,
        &mut contenders,
        successors_inclusive,
        predecessors,
        active,
        &get_data,
        functional,
        cache,
    );

    let chosen = contenders.iter().max_by_key(|set| set.len()).unwrap();
    (chosen.clone(), best_value)
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn dfs<
    D: Direction,
    const B: usize,
    F: CauchyMeanValueFunctional,
    G: Fn(usize) -> I + Clone,
    I: Iterator<Item = Observation<(), F::Response, F::Censoring>> + Clone,
>(
    in_set: BitSet<B>,
    forbidden: &BitSet<B>,
    best_value: &mut f64,
    contenders: &mut Vec<BitSet<B>>,
    successors_inclusive: &[BitSet<B>],
    predecessors: &[BitSet<B>],
    active: &BitSet<B>,
    get_data: &G,
    functional: &F,
    cache: &mut HashMap<BitSet<B>, f64>,
) {
    let decided = in_set.union(forbidden);
    let undecided = active.difference(&decided);

    // Available v: undecided and all predecessors (incl) of v are already in in_set
    // (equivalently, strict preds ⊆ in_set since v itself is not in_set yet).
    match undecided.min() {
        None => {
            // We found a lower set that we should check
            let lower_set = in_set;
            if lower_set.is_empty() {
                return;
            }

            let value = functional.evaluate_partial_order(
                &lower_set,
                successors_inclusive,
                predecessors,
                get_data,
                cache,
            );
            let is_better = if D::IS_INCREASING {
                value < *best_value
            } else {
                value > *best_value
            };
            if is_better {
                *best_value = value;
                contenders.clear();
                contenders.push(lower_set);
            } else if value == *best_value {
                contenders.retain(|existing| !existing.is_superset_of(&lower_set));
                contenders.push(lower_set);
            }
        }
        Some(minimal_element) => {
            // include v
            dfs::<D, _, _, _, _>(
                in_set.with_singleton(minimal_element),
                forbidden,
                best_value,
                contenders,
                successors_inclusive,
                predecessors,
                active,
                get_data,
                functional,
                cache,
            );

            // exclude v => forbid v and everything above it (keeps "lower set" property)
            dfs::<D, _, _, _, _>(
                in_set,
                &forbidden.union(&successors_inclusive[minimal_element]),
                best_value,
                contenders,
                successors_inclusive,
                predecessors,
                active,
                get_data,
                functional,
                cache,
            );
        }
    }
}

#[cfg(test)]
mod test {
    use crate::functionals::Average;
    use crate::partial_order::functionals::algorithm_pre_sorted;
    use crate::structures::{Increasing, Observation};

    #[test]
    fn chain() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [1.0, 3.0, 2.0].into_iter(),
            &[(0, 1), (1, 2)],
            &Average,
        );
        let expected = vec![1.0, 2.5, 2.5];
        assert_eq!(result, expected);
    }

    #[test]
    fn chain_antitone() {
        // 0 <= 1 <= 2 <= 3, strictly decreasing, all pool to the mean 2.5
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [4.0, 3.0, 2.0, 1.0].into_iter(),
            &[(0, 1), (1, 2), (2, 3)],
            &Average,
        );
        let expected = vec![2.5; 4];
        assert_eq!(result, expected);
    }

    #[test]
    fn chain_semitone() {
        // 0 <= 1 <= 2 <= 3, only a local violation at (1, 2)
        // PAV on a chain gives [1, 3, 3, 3]
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [1.0, 4.0, 2.0, 3.0].into_iter(),
            &[(0, 1), (1, 2), (2, 3)],
            &Average,
        );
        let expected = vec![1.0, 3.0, 3.0, 3.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn split() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [1.0, 2.0, 0.0].into_iter(),
            &[(0, 1), (0, 2)],
            &Average,
        );
        let expected = vec![0.5, 2.0, 0.5];
        assert_eq!(result, expected);
    }

    #[test]
    fn split_asymmetric() {
        assert_eq!(
            algorithm_pre_sorted::<Increasing, _, _>(
                [6.0, 12.0, 12.0, 24.0].into_iter(),
                &[(3, 2), (3, 1), (1, 0)],
                &Average,
            ),
            vec![13.5; 4],
        );
    }

    #[test]
    fn join() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [2.0, 3.0, 1.0].into_iter(),
            &[(1, 0), (1, 2)],
            &Average,
        );
        assert_eq!(result, vec![2.0, 2.0, 2.0]);
    }

    #[test]
    fn join_unequal() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [20.0, 10.0, 0.0].into_iter(),
            &[(0, 2), (1, 2)],
            &Average,
        );
        assert_eq!(result, vec![10.0; 3]);
    }

    #[test]
    fn diamond_mixed() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [3.0, 1.0, 4.0, 2.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 3)],
            &Average,
        );
        let expected = vec![2.0, 2.0, 3.0, 3.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn diamond_antitone() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [5.0, 1.0, 1.0, 1.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 3)],
            &Average,
        );
        let expected = vec![2.0; 4];
        assert_eq!(result, expected);
    }

    #[test]
    fn diamond_semitone() {
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [4.0, 6.0, 2.0, 3.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 3)],
            &Average,
        );
        let expected = vec![3.0, 4.5, 3.0, 4.5];
        assert_eq!(result, expected);
    }

    #[test]
    fn branched_isotone() {
        // A small branched poset:
        //    3
        //   / \
        //  1   2   -> 4
        //   \ /
        //    0
        //
        // Edges: 0 <= 1, 0 <= 2, 1 <= 3, 2 <= 3, 2 <= 4
        // Values already isotone.
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [1.0, 2.0, 2.5, 3.0, 4.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 3), (2, 4)],
            &Average,
        );
        let expected = vec![1.0, 2.0, 2.5, 3.0, 4.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn branched_semitone() {
        // A small branched poset:
        //    3
        //   / \
        //  1   2   -> 4
        //   \ /
        //    0
        //
        // Edges: 0 <= 1, 0 <= 2, 1 <= 3, 2 <= 3, 2 <= 4
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [1.0, 2.0, 3.0, 4.0, 0.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 3), (2, 4)],
            &Average,
        );
        let expected = vec![1.0, 2.0, 1.5, 4.0, 1.5];
        assert_eq!(result, expected);
    }

    #[test]
    fn branched_semitone_merge_further() {
        // A small branched poset:
        //    3
        //   / \
        //  1   2   -> 4
        //   \ /
        //    0
        //
        // Edges: 0 <= 1, 0 <= 2, 1 <= 3, 2 <= 3, 2 <= 4
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [3.0, 6.0, 3.0, 12.0, 0.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 3), (2, 4)],
            &Average,
        );
        let expected = vec![2.0, 6.0, 2.0, 12.0, 2.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn branched_complicated() {
        // A "two-chain ladder" poset where a naive
        // linear-extension PAV is suboptimal:
        //
        //        5
        //       / \
        //      3   4
        //     /     \
        //    1       2
        //     \     /
        //        0
        //
        // Edges: 0 <= 1, 0 <= 2, 1 <= 3, 2 <= 4, 3 <= 5, 4 <= 5
        //
        // If you take a topological order like [0, 1, 2, 3, 4, 5]
        // and just run 1D PAV on that order, you'll get something like
        //   [0.0, 0.0, 20.0/3.0, 20.0/3.0, 20.0/3.0, 10.0]
        // which *is* isotone on that total order but *not* optimal
        // for the partial order.
        //
        // The correct max–min / isotone solution on the poset is
        //   [0.0, 0.0, 5.0, 10.0, 5.0, 10.0]
        // where only nodes 2 and 4 are pooled.
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [0.0, 0.0, 10.0, 10.0, 0.0, 10.0].into_iter(),
            &[(0, 1), (0, 2), (1, 3), (2, 4), (3, 5), (4, 5)],
            &Average,
        );
        let expected = vec![0.0, 0.0, 5.0, 10.0, 5.0, 10.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn three_node_star_gpav_counterexample() {
        // Star poset: 0 <= 1, 0 <= 2
        // responses y = [8, 7, 0], unit weights
        //
        // True isotonic regression (via the max–min characterization
        // or KKT conditions) is [4, 7, 4] with SSE = 32.
        // The current GPAV-style implementation, with edges ordered
        // as below, topologically sorts as [0, 1, 2] and returns
        // [5, 5, 5] with SSE = 38, which is suboptimal.
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [8.0, 7.0, 0.0].into_iter(),
            &[(0, 1), (0, 2)],
            &Average,
        );
        let expected = vec![4.0, 7.0, 4.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn disconnected() {
        let result =
            algorithm_pre_sorted::<Increasing, _, _>([3.0, 2.0, 1.0].into_iter(), &[], &Average);
        let expected = vec![3.0, 2.0, 1.0];
        assert_eq!(result, expected);
    }

    #[test]
    fn two_node_weighted_antitone() {
        // 0 <= 1, y = [10, 0], w = [10, 1]
        // Violates isotonicity, so both pool to weighted mean:
        // t = (10 * 10 + 1 * 0) / (10 + 1) = 100 / 11
        let result = algorithm_pre_sorted::<Increasing, _, _>(
            [
                Observation {
                    x: (),
                    y: 10.0,
                    observed: (),
                    weight: 10.0,
                },
                Observation {
                    x: (),
                    y: 0.0,
                    observed: (),
                    weight: 1.0,
                },
            ]
            .into_iter(),
            &[(0, 1)],
            &Average,
        );
        let expected = vec![100.0 / 11.0; 2];
        assert_eq!(result, expected);
    }
}
