use crate::partial_order::structures::{NodeID, OrderingInfo};
use crate::partial_order::{BitSet, PredictionWorkspace};
use itertools::Itertools;
use std::collections::VecDeque;
use std::collections::hash_set::Drain;

#[must_use]
pub fn derive_transitive_reduction(
    covariates: &[f64], // point-major: x_i is covariates[i*covariate_dimension .. (i+1)*covariate_dimension]
    n_covariate: usize,
    covariate_dimension: usize,
) -> Vec<(usize, usize)> {
    debug_assert_ne!(n_covariate, 0);
    debug_assert_ne!(covariate_dimension, 0);
    debug_assert_eq!(covariates.len(), n_covariate * covariate_dimension);
    debug_assert!(
        covariates
            .chunks_exact(covariate_dimension)
            .tuple_windows()
            .all(|(l, r)| {
                for (ll, rr) in l.iter().zip(r.iter()) {
                    if ll < rr {
                        return true;
                    }
                    if ll > rr {
                        return false;
                    }
                }
                // not unique
                false
            }),
        "Covariates must be lexically sorted and deduplicated",
    );

    // Helper to view the i-th point as a slice.
    #[inline]
    fn row(covariates: &[f64], i: usize, dim: usize) -> &[f64] {
        &covariates[i * dim..(i + 1) * dim]
    }

    // Strict dominance test: x_u < x_v (coordinatewise <= and at least one strict).
    let dominates_strict = |u: usize, v: usize| -> bool {
        let (xu, xv) = (
            row(covariates, u, covariate_dimension),
            row(covariates, v, covariate_dimension),
        );
        let mut strictly_less = false;
        for (&a, &b) in xu.iter().zip(xv.iter()) {
            if a > b {
                return false; // violates dominance in this coordinate
            }
            if a < b {
                strictly_less = true;
            }
        }
        strictly_less
    };

    // 2) Bitset reachability. R[i] is the bitset of vertices reachable from vertex i (by input index).
    let n = n_covariate;
    let words = n.div_ceil(64);
    let mut reach = vec![vec![0_u64; words]; n];

    // Bit ops
    #[inline]
    fn test_bit(bits: &[u64], idx: usize) -> bool {
        let w = idx >> 6;
        let b = idx & 63;
        (bits[w] >> b) & 1u64 == 1
    }
    #[inline]
    fn set_bit(bits: &mut [u64], idx: usize) {
        let w = idx >> 6;
        let b = idx & 63;
        bits[w] |= 1u64 << b;
    }
    #[inline]
    fn or_inplace(dst: &mut [u64], src: &[u64]) {
        for (d, s) in dst.iter_mut().zip(src.iter()) {
            *d |= *s;
        }
    }

    // Final cover edges use original indices.
    let mut edges: Vec<(usize, usize)> = Vec::with_capacity(n);

    // 3) Reverse topological pass over the identity order (already lex-sorted).
    // For each u from n-1 down to 0, scan v in (u+1..n).
    //   - If u ≺ v and v already in R_acc => redundant
    //   - Else keep edge (u,v), accumulate reach[v] into R_acc, and set v in R_acc
    for upos in (0..n).rev() {
        // Mutably borrow reach[upos] and immutably access the tail via split_at_mut to satisfy the borrow checker.
        let (prefix, suffix) = reach.split_at_mut(upos + 1);
        let r_acc = &mut prefix[upos];
        // zero the accumulator for this u
        for w in r_acc.iter_mut() {
            *w = 0;
        }

        for vpos in (upos + 1)..n {
            if dominates_strict(upos, vpos) {
                if test_bit(r_acc, vpos) {
                    continue; // redundant
                }
                edges.push((upos, vpos));
                // reach[vpos] lives in the suffix slice at index vpos - (upos + 1)
                let rv = &suffix[vpos - (upos + 1)];
                or_inplace(r_acc, rv);
                set_bit(r_acc, vpos);
            }
        }
        // r_acc already stored in reach[upos]
    }

    edges
}

/// Compute the lower or upper neighbors of a new target.
///
/// The comments are written for when computing the lower neighbors.
pub fn find_neighbors<'a>(
    target: &[f64],
    covariates: &[f64],
    ordering_info: &OrderingInfo,
    workspace: &'a mut PredictionWorkspace,
    lower: bool,
) -> Drain<'a, NodeID> {
    debug_assert!(workspace.stack.is_empty());
    debug_assert!(workspace.visited.is_empty());
    debug_assert!(workspace.found.is_empty());

    let covariate_dimension = target.len();
    let get_covariate =
        |i: NodeID| &covariates[i * covariate_dimension..(i + 1) * covariate_dimension];

    type Comparison = fn(&[f64], &[f64]) -> bool;
    let (comparison, start_set, next) = if lower {
        (
            point_ge_point as Comparison,
            &ordering_info.max,
            &ordering_info.smaller,
        )
    } else {
        (
            point_le_point as Comparison,
            &ordering_info.min,
            &ordering_info.larger,
        )
    };

    workspace.stack.extend_from_slice(start_set);

    while let Some(index) = workspace.stack.pop() {
        if workspace.visited.contains(&index) {
            continue;
        }
        workspace.visited.insert(index);

        if comparison(target, get_covariate(index)) {
            // We found a point that is fully below our target, so stop the recursion

            // Remove any previously-found children that are smaller than this one
            // TODO: Costs O(capacity), can we be faster?
            workspace
                .found
                .retain(|&child| !comparison(get_covariate(index), get_covariate(child)));
            workspace.found.insert(index);
        } else {
            // Our index is still above our target or the points are incomparable, continue downward
            // TODO: Try extending with the entire slice, filtering happens after anyway
            workspace.stack.extend(
                next.get(index)
                    .iter()
                    .filter(|&&i| !workspace.visited.contains(&i)),
            );
        }
    }
    debug_assert!(workspace.stack.is_empty());
    workspace.visited.clear();
    workspace.found.drain()
}

#[inline]
fn point_ge_point(a: &[f64], b: &[f64]) -> bool {
    // a >= b component-wise (non-strict)
    a.iter().zip(b.iter()).all(|(x, y)| x >= y)
}

#[inline]
fn point_le_point(a: &[f64], b: &[f64]) -> bool {
    // a <= b component-wise (non-strict)
    a.iter().zip(b.iter()).all(|(x, y)| x <= y)
}

#[cfg(test)]
mod tests {
    use super::derive_transitive_reduction;
    use crate::routines::{argsort_unstable_by, lexicographic_cmp};
    use crate::structures::Increasing;
    use std::collections::HashSet;

    fn lexically_sort(covariates: &[f64], n: usize, d: usize) -> Vec<f64> {
        assert_eq!(covariates.len(), n * d);

        let get_cdf = |i| &covariates[i * d..(i + 1) * d];
        let order = argsort_unstable_by::<Increasing, _>(
            |a, b| lexicographic_cmp(get_cdf(a), get_cdf(b)),
            n,
        );

        order
            .into_iter()
            .map(|i| get_cdf(i).iter())
            .flatten()
            .copied()
            .collect()
    }

    fn edges_sorted(mut e: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        e.sort();
        e
    }

    #[test]
    fn small_2d_example() {
        // Points:
        // 0: (0,0)
        // 1: (0,1)
        // 2: (1,0)
        // 3: (1,1)
        // Covers: (0,0) is covered by (0,1) and (1,0); (0,1) and (1,0) are covered by (1,1)
        // Edges: (0->1), (0->2), (1->3), (2->3)
        let pts = [
            0.0, 0.0, // 0
            0.0, 1.0, // 1
            1.0, 0.0, // 2
            1.0, 1.0, // 3
        ];
        let edges = derive_transitive_reduction(&pts, 4, 2);
        assert_eq!(
            edges_sorted(edges),
            edges_sorted(vec![(0, 1), (0, 2), (1, 3), (2, 3)])
        );
    }

    #[test]
    fn mixed_3d_example() {
        // 0: (0,0,0)
        // 1: (0,1,1)
        // 2: (1,0,1)
        // 3: (1,1,0)
        // 4: (1,1,1)
        // Covers: 0 is covered by 1,2,3; and 1,2,3 are covered by 4
        let pts = [
            0.0, 0.0, 0.0, // 0
            0.0, 1.0, 1.0, // 1
            1.0, 0.0, 1.0, // 2
            1.0, 1.0, 0.0, // 3
            1.0, 1.0, 1.0, // 4
        ];
        let edges = derive_transitive_reduction(&pts, 5, 3);
        assert_eq!(
            edges_sorted(edges),
            edges_sorted(vec![(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4)])
        );
    }

    // ------------------------
    // Small utilities for tests
    // ------------------------

    fn row<'a>(covs: &'a [f64], d: usize, i: usize) -> &'a [f64] {
        &covs[i * d..i * d + d]
    }

    fn sorted_pairs(mut v: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        v.sort_unstable();
        v
    }

    // Naive oracle: strict component-wise order and cover edges in O(n^3)
    fn leq(a: &[f64], b: &[f64]) -> bool {
        a.iter().zip(b).all(|(x, y)| x <= y)
    }
    fn lt(a: &[f64], b: &[f64]) -> bool {
        leq(a, b) && a.iter().zip(b).any(|(x, y)| x < y)
    }
    fn naive_covers(covs: &[f64], n: usize, d: usize) -> Vec<(usize, usize)> {
        let mut e = Vec::new();
        for i in 0..n {
            let ai = row(covs, d, i);
            for j in 0..n {
                if i == j {
                    continue;
                }
                let aj = row(covs, d, j);
                if lt(ai, aj) {
                    let mut covered = false;
                    for k in 0..n {
                        if k == i || k == j {
                            continue;
                        }
                        let ak = row(covs, d, k);
                        if lt(ai, ak) && lt(ak, aj) {
                            covered = true;
                            break;
                        }
                    }
                    if !covered {
                        e.push((i, j));
                    }
                }
            }
        }
        sorted_pairs(e)
    }

    #[test]
    fn derive_constraints_grid_example() {
        // Same as minimal_dominators_small_posets: expect cover edges
        let covs = vec![
            0.0, 0.0, // 0
            0.0, 1.0, // 1
            1.0, 0.0, // 2
            1.0, 1.0, // 3
            2.0, 2.0, // 4
        ];
        let edges = sorted_pairs(derive_transitive_reduction(&covs, 5, 2));
        let expected = sorted_pairs(vec![(0, 1), (0, 2), (1, 3), (2, 3), (3, 4)]);
        assert_eq!(edges, expected);
    }

    #[test]
    fn derive_constraints_chain_and_antichain() {
        // Chain in 3D: 0<1<2<3<4
        let mut covs = Vec::new();
        for i in 0..5 {
            covs.extend_from_slice(&[i as f64, (i * 2) as f64, (i * 3) as f64]);
        }
        let edges_chain = sorted_pairs(derive_transitive_reduction(&covs, 5, 3));
        let expected_chain = sorted_pairs(vec![(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert_eq!(edges_chain, expected_chain);

        // Antichain in 2D: pairwise incomparable
        let covs2 = vec![0.0, 1.0, 0.2, 0.9, 0.5, 0.5, 1.0, 0.0];
        let edges_anti = derive_transitive_reduction(&covs2, 4, 2);
        assert!(edges_anti.is_empty());
    }

    #[test]
    fn cross_check_against_naive_on_random_small() {
        // Deterministic pseudo-random generator (no rand crate)
        fn lcg(seed: &mut u64) -> u32 {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (*seed >> 32) as u32
        }
        fn next_f(seed: &mut u64) -> f64 {
            let v = lcg(seed) as f64 / (u32::MAX as f64);
            // keep values on a small grid to create comparable ties sometimes
            (v * 10.0).floor() / 10.0
        }
        let mut seed = 123456789u64;

        for _ in 0..5 {
            let n = 10usize;
            let d = 3usize;
            let mut covs = Vec::with_capacity(n * d);
            for _ in 0..n {
                for _ in 0..d {
                    covs.push(next_f(&mut seed));
                }
            }
            let covs = lexically_sort(&covs, n, d);
            let got = sorted_pairs(derive_transitive_reduction(&covs, n, d));
            let want = naive_covers(&covs, n, d);
            assert_eq!(got, want, "Mismatch vs naive cover graph");
        }
    }

    #[test]
    fn small_2d() {
        // Points:
        // 0: (0,0)
        // 1: (0,1)
        // 2: (1,0)
        // 3: (1,1)
        // 4: (2,2)
        let covariates = vec![0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 2.0, 2.0];
        let mut edges = derive_transitive_reduction(&covariates, 5, 2);

        // Expected covers: 0->1, 0->2, 1->3, 2->3, 3->4
        edges.sort_unstable();
        let expected = vec![(0, 1), (0, 2), (1, 3), (2, 3), (3, 4)];
        assert_eq!(edges, expected);
    }

    // -----------------------------
    // Helpers
    // -----------------------------

    // Verify generic invariants that should hold for any valid output.
    fn assert_constraint_invariants(
        covariates: &[f64],
        n_covariate: usize,
        cov_dim: usize,
        constraints: &[(usize, usize)],
    ) {
        // 1) Indices must be in range, lo < hi, and unique.
        let mut seen = HashSet::new();
        for &(lo, hi) in constraints {
            assert!(
                lo < n_covariate && hi < n_covariate,
                "Indices must be < n_covariate. Got ({lo}, {hi}) with n_covariate={n_covariate}."
            );
            assert!(
                seen.insert((lo, hi)),
                "Duplicate constraint ({lo}, {hi}) found; constraints should be unique."
            );

            // 2) Soundness: every returned pair must be element-wise <=.
            let lo_row = row(covariates, cov_dim, lo);
            let hi_row = row(covariates, cov_dim, hi);
            assert!(
                leq(lo_row, hi_row),
                "Unsound constraint: row {lo:?} = {lo_row:?} is not <= row {hi:?} = {hi_row:?} element-wise."
            );
        }
    }

    // Build a covariate-major flattened buffer from rows.
    fn flatten<const N: usize, const M: usize>(rows: &[&[f64; N]; M]) -> Vec<f64> {
        rows.iter().flat_map(|r| r.iter().copied()).collect()
    }

    // Is (lo, hi) in the returned constraints?
    fn contains_pair(pairs: &[(usize, usize)], lo: usize, hi: usize) -> bool {
        pairs.iter().any(|&(i, j)| i == lo && j == hi)
    }

    #[test]
    fn one_covariate_any_dim() {
        let rows = &[&[1.0, -2.0, 3.0]];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 1, 3);
        assert!(
            out.is_empty(),
            "With a single covariate, there cannot be any pairwise constraints."
        );
    }

    #[test]
    fn one_dimensional_values_soundness() {
        // d = 1 reduces to ordinary numeric ordering.
        let rows = &[&[1.0], &[2.0], &[3.0], &[5.0]];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 4, 1);
        assert_constraint_invariants(&cov, 4, 1, &out);
    }

    // -----------------------------
    // Structural correctness on crafted data
    // -----------------------------

    #[test]
    fn strict_chain_requires_consecutive_edges() {
        // Rows form a strict chain: r0 <= r1 <= r2 <= r3 <= r4 element-wise.
        // Any correct algorithm (cover edges or all pairs) must include (i, i+1).
        let rows = &[
            &[0.0, 0.0, 0.0],
            &[1.0, 1.0, 1.0],
            &[2.0, 2.0, 2.0],
            &[3.0, 3.0, 3.0],
            &[4.0, 4.0, 4.0],
        ];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 5, 3);
        assert_constraint_invariants(&cov, 5, 3, &out);

        for i in 0..4 {
            assert!(
                contains_pair(&out, i, i + 1),
                "Chain must contain consecutive edge ({}, {}). Got: {:?}",
                i,
                i + 1,
                out
            );
        }
        // Upper bound sanity check: cannot exceed all-pairs.
        assert!(out.len() <= 5 * 4 / 2);
    }

    #[test]
    fn anti_chain_yields_no_constraints() {
        // No two rows are mutually comparable element-wise (strict antichain).
        let rows = &[&[0.0, 2.0], &[1.0, 1.0], &[2.0, 0.0]];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 3, 2);
        // Soundness + expect empty: any pair would violate element-wise <=
        assert!(
            out.is_empty(),
            "Antichain should yield no constraints, got {out:?}"
        );
    }

    #[test]
    fn indexing_is_row_major_covariate_major() {
        // Construct rows so that only row-major interpretation forms a chain.
        // If an implementation mistakenly indexes as column-major, the chain property would likely break.
        let rows = &[
            // r0 <= r1 <= r2 (element-wise)
            &[1.0, 10.0, 100.0],
            &[2.0, 20.0, 200.0],
            &[3.0, 30.0, 300.0],
        ];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 3, 3);
        assert_constraint_invariants(&cov, 3, 3, &out);
        // Must contain consecutive edges if the rows were read correctly.
        assert!(
            contains_pair(&out, 0, 1),
            "Missing (0,1) suggests wrong flattening/stride."
        );
        assert!(
            contains_pair(&out, 1, 2),
            "Missing (1,2) suggests wrong flattening/stride."
        );
    }

    #[test]
    fn handles_negative_and_large_values() {
        let rows = &[&[-1e12, -2.0, 0.0], &[0.0, 0.0, 0.0], &[1e12, 2.0, 1e3]];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 3, 3);
        assert_constraint_invariants(&cov, 3, 3, &out);
    }

    #[test]
    fn handles_infinities() {
        let rows = &[
            &[0.0, 0.0],
            &[f64::INFINITY, 0.0],
            &[f64::INFINITY, f64::INFINITY],
        ];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 3, 2);
        assert_constraint_invariants(&cov, 3, 2, &out);
        // Should contain at least the chain edges (0,1) and (1,2)
        assert!(contains_pair(&out, 0, 1));
        assert!(contains_pair(&out, 1, 2));
    }

    #[test]
    fn medium_stress_no_panic() {
        // Not about exact content; just ensure it scales and respects invariants for moderately large inputs.
        let n = 2000;
        let d = 3;
        let mut cov = Vec::with_capacity(n * d);
        for i in 0..n {
            // Increasing chain so work is predictable
            for _ in 0..d {
                cov.push(i as f64);
            }
        }
        let out = derive_transitive_reduction(&cov, n, d);
        assert_constraint_invariants(&cov, n, d, &out);
    }

    // a <= b element-wise with at least one strict <
    fn strict_leq(a: &[f64], b: &[f64]) -> bool {
        let mut any_strict = false;
        for (x, y) in a.iter().zip(b.iter()) {
            if x > y {
                return false;
            }
            if x < y {
                any_strict = true;
            }
        }
        any_strict
    }

    // Is there an element k strictly between i and j?
    fn has_middle(cov: &[f64], n: usize, d: usize, i: usize, j: usize) -> bool {
        let ri = row(cov, d, i);
        let rj = row(cov, d, j);
        (0..n).any(|k| {
            if k == i || k == j {
                return false;
            }
            let rk = row(cov, d, k);
            strict_leq(ri, rk) && strict_leq(rk, rj)
        })
    }

    // Assert every returned pair is a cover (and properly ordered as i<j and i<j in poset)
    fn assert_only_cover_edges(cov: &[f64], n: usize, d: usize, edges: &[(usize, usize)]) {
        for &(i, j) in edges {
            assert!(i < j, "Pairs must be normalized with i < j, got ({i},{j}).");
            assert!(i < n && j < n, "Index out of bounds ({i},{j}) for n={n}.");

            let ri = row(cov, d, i);
            let rj = row(cov, d, j);

            assert!(
                strict_leq(ri, rj),
                "Returned pair ({i},{j}) is not ordered in the poset."
            );
            assert!(
                !has_middle(cov, n, d, i, j),
                "Returned pair ({i},{j}) is transitive (there exists k with i<k<j); not a cover."
            );
        }
    }

    fn contains(edges: &[(usize, usize)], i: usize, j: usize) -> bool {
        edges.iter().any(|&(u, v)| u == i && v == j)
    }

    // 1) Strict chain: only consecutive edges are covers. Any (i,j) with j - i > 1 is not a cover.
    #[test]
    fn chain_covers_only() {
        // r0 < r1 < r2 < r3 < r4
        let rows = &[
            &[0.0, 0.0],
            &[1.0, 1.0],
            &[2.0, 2.0],
            &[3.0, 3.0],
            &[4.0, 4.0],
        ];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 5, 2);

        // Must be only covers
        assert_only_cover_edges(&cov, 5, 2, &out);

        // In a chain, covers are exactly consecutive pairs
        for i in 0..4 {
            assert!(
                contains(&out, i, i + 1),
                "Missing cover edge ({},{})",
                i,
                i + 1
            );
        }
        for i in 0..5 {
            for j in i + 2..5 {
                assert!(
                    !contains(&out, i, j),
                    "Found non-cover transitive edge ({},{}) in a chain.",
                    i,
                    j
                );
            }
        }
    }

    // 2) Diamond poset: 0 < 1, 0 < 2, 1 < 3, 2 < 3 are covers; 0 < 3 is transitive and must be absent.
    #[test]
    fn diamond_covers_only() {
        // 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
        let rows = &[
            &[0.0, 0.0], // 0
            &[0.0, 1.0], // 1
            &[1.0, 0.0], // 2
            &[1.0, 1.0], // 3
        ];
        let cov = flatten(rows);
        let out = derive_transitive_reduction(&cov, 4, 2);

        assert_only_cover_edges(&cov, 4, 2, &out);

        // Expected covers
        for &(i, j) in &[(0, 1), (0, 2), (1, 3), (2, 3)] {
            assert!(contains(&out, i, j), "Missing cover edge ({},{})", i, j);
        }
        // Must not contain the long edge 0->3
        assert!(
            !contains(&out, 0, 3),
            "Found transitive (non-cover) edge (0,3) in a diamond."
        );
        // 1 and 2 are incomparable
        assert!(
            !contains(&out, 1, 2) && !contains(&out, 2, 1),
            "Incomparable nodes (1,2) should not have an edge."
        );
    }

    // 3) Boolean lattice B3: covers increase Hamming weight by exactly 1.
    // Assert the output set is a subset of the true cover set; optionally assert equality.
    #[test]
    fn boolean_lattice_b3_only_covers() {
        // Build all 3-bit 0/1 vectors, ordered by index 0..7
        let mut rows: Vec<Vec<f64>> = Vec::new();
        for x in 0u8..8 {
            let v = [
                ((x >> 0) & 1) as f64,
                ((x >> 1) & 1) as f64,
                ((x >> 2) & 1) as f64,
            ];
            rows.push(v.to_vec());
        }
        let cov: Vec<f64> = rows.iter().flatten().copied().collect();
        let n = 8;
        let d = 3;

        // Compute expected cover set for Boolean lattice: A ⊂ B and |B|-|A| = 1
        let mut expected: HashSet<(usize, usize)> = HashSet::new();
        let weight = |x: u8| x.count_ones() as i32;
        for a in 0u8..8 {
            for b in 0u8..8 {
                if a == b {
                    continue;
                }
                if (a & b) == a && weight(b) - weight(a) == 1 {
                    expected.insert((a as usize, b as usize));
                }
            }
        }

        let cov = lexically_sort(&cov, n, d);
        let out = derive_transitive_reduction(&cov, n, d);
        assert_only_cover_edges(&cov, n, d, &out);

        // Ensure every returned edge is indeed a cover in the lattice
        for &(i, j) in &out {
            assert!(
                expected.contains(&(i, j)),
                "Returned edge ({i},{j}) is not a cover in B3."
            );
        }

        // If your spec requires completeness (return all covers), also assert equality:
        // let out_set: HashSet<_> = out.iter().copied().collect();
        // assert_eq!(out_set, expected, "Output must be exactly the cover set of B3.");
    }
}

/// Compute strict successor closure from Hasse edges (DAG).
pub fn compute_transitive_closure<const B: usize>(
    n: usize,
    edges: &[(usize, usize)],
) -> (Vec<BitSet<B>>, Vec<BitSet<B>>, Vec<usize>, Vec<usize>) {
    assert!(n <= BitSet::<B>::capacity());

    // Build adjacency + indegree for topo sort
    let mut out = vec![Vec::<usize>::new(); n];
    let mut in_degree = vec![0usize; n];

    for &(i, j) in edges {
        assert!(i < n && j < n);
        out[i].push(j);
        in_degree[j] += 1;
    }

    // Kahn topo sort
    let mut queue = VecDeque::new();
    for (v, &degree) in in_degree.iter().enumerate() {
        if degree == 0 {
            queue.push_back(v);
        }
    }
    let mut topological_order = Vec::with_capacity(n);
    while let Some(v) = queue.pop_front() {
        topological_order.push(v);
        for &w in &out[v] {
            in_degree[w] -= 1;
            if in_degree[w] == 0 {
                queue.push_back(w);
            }
        }
    }
    // If this is truly a Hasse diagram, it must be acyclic.
    assert_eq!(
        topological_order.len(),
        n,
        "graph has a cycle; not a Hasse diagram / not a DAG",
    );
    debug_assert_eq!(topological_order.iter().sum::<usize>(), (n - 1) * n / 2);
    debug_assert!(in_degree.into_iter().all(|d| d == 0));

    let reverse_topological_order = {
        let mut x = vec![0; n];
        for (i, &u) in topological_order.iter().enumerate() {
            x[u] = i;
        }
        debug_assert_eq!(x.iter().sum::<usize>(), (n - 1) * n / 2);
        x
    };

    // DP in reverse topo: succ[v] |= {w} ∪ succ[w] for each edge v->w
    let mut successors = vec![BitSet::new(); n];
    let mut predecessors = vec![BitSet::new(); n];
    for topological_u in (0..n).rev() {
        let u = topological_order[topological_u];
        for &v in &out[u] {
            let topological_v = reverse_topological_order[v];
            let tmp = successors[topological_u].clone();
            successors[topological_u] = successors[topological_u]
                .with_singleton(topological_v)
                .union(&successors[topological_v]);

            let delta = successors[topological_u].difference(&tmp);
            // only process newly added reachables
            for topological_successor in delta.iter() {
                predecessors[topological_successor].add(topological_u);
            }
        }
    }

    (
        successors,
        predecessors,
        topological_order,
        reverse_topological_order,
    )
}
