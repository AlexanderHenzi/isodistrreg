use crate::structures::Direction;
use osqp::CscMatrix;
use std::borrow::Cow;

#[allow(clippy::doc_overindented_list_items)]
/// Update A in-place to encode HRO constraints per edge-row.
/// For each row r that originally encoded x_i <= x_j:
/// - If D::IS_INCREASING:
///     HRO:  (x_i / S_i) - (x_j / S_j) >= 0, if S_i > 0 and S_j > 0
///     else:            x_j - x_i       >= 0
/// - Else (decreasing):
///     HRO:  (x_j / S_j) - (x_i / S_i) >= 0, if S_i > 0 and S_j > 0
///     else:            x_i - x_j       >= 0
///
/// Invariants expected:
/// - existing_matrix has exactly two nonzeros per row (built by your builder or
///   previously modified only by this function).
/// - survival.len() == existing_matrix.ncols.
/// - Within each row r, the two variables (i, j) satisfy:
///     if D::IS_INCREASING      then survival[i] <= survival[j]
///     if !D::IS_INCREASING     then survival[i] >= survival[j]
pub fn update_constraint_matrix<D: Direction>(existing_matrix: &mut CscMatrix, survival: &[f64]) {
    let m = existing_matrix.nrows;
    let n = existing_matrix.ncols;

    assert_eq!(
        survival.len(),
        n,
        "survival length ({}) must equal number of columns ({})",
        survival.len(),
        n
    );

    // Borrow CSC slices. We'll only mutate data.
    let indptr = match &existing_matrix.indptr {
        Cow::Borrowed(s) => *s,
        Cow::Owned(v) => v.as_slice(),
    };
    let indices = match &existing_matrix.indices {
        Cow::Borrowed(s) => *s,
        Cow::Owned(v) => v.as_slice(),
    };
    let data = existing_matrix.data.to_mut();

    // Temporary storage: remember first (col, position) seen for each row.
    // Using usize::MAX as "none" sentinel avoids Option overhead.
    let mut first_col = vec![usize::MAX; m];
    let mut first_pos = vec![0usize; m];

    // Helper: choose (i, j) from two columns using the survival ordering guarantee.
    #[inline]
    fn orient_pair<D: Direction>(
        ca: usize,
        sa: f64,
        cb: usize,
        sb: f64,
    ) -> (usize, f64, usize, f64) {
        if D::IS_INCREASING {
            // i has smaller (or equal) survival, j the larger.
            // Tie-break deterministically by column index.
            if sa < sb || (sa == sb && ca <= cb) {
                (ca, sa, cb, sb) // i = a, j = b
            } else {
                (cb, sb, ca, sa) // i = b, j = a
            }
        } else {
            // Decreasing: i has larger (or equal) survival, j the smaller.
            if sa > sb || (sa == sb && ca >= cb) {
                (ca, sa, cb, sb) // i = a, j = b
            } else {
                (cb, sb, ca, sa) // i = b, j = a
            }
        }
    }

    // Pass over columns; when we see the second entry for a row, we have the full pair.
    for col in 0..n {
        #[allow(clippy::needless_range_loop)]
        for pos in indptr[col]..indptr[col + 1] {
            let r = indices[pos];
            debug_assert!(r < m, "row index out of range");

            if first_col[r] == usize::MAX {
                first_col[r] = col;
                first_pos[r] = pos;
                continue;
            }

            // We found the second column for this row r.
            let col_a = first_col[r];
            let pos_a = first_pos[r];
            let col_b = col;
            let pos_b = pos;

            // Identify i and j using the survival guarantee.
            let s_a = survival[col_a];
            let s_b = survival[col_b];
            let (i_col, s_i, _j_col, s_j) = orient_pair::<D>(col_a, s_a, col_b, s_b);

            // Map columns back to their storage positions for this row.
            let (pos_i, pos_j) = if i_col == col_a {
                (pos_a, pos_b)
            } else {
                (pos_b, pos_a)
            };

            // Decide HRO vs fallback (simple order) for this row.
            let use_hro = s_i > 0.0 && s_j > 0.0;

            if use_hro {
                // HRO coefficients: signs depend on direction.
                if D::IS_INCREASING {
                    // (x_i/S_i) - (x_j/S_j) >= 0
                    data[pos_i] = 1.0 / s_i;
                    data[pos_j] = -1.0 / s_j;
                } else {
                    // (x_j/S_j) - (x_i/S_i) >= 0
                    data[pos_i] = -1.0 / s_i;
                    data[pos_j] = 1.0 / s_j;
                }
            } else {
                // Fallback to simple order.
                if D::IS_INCREASING {
                    // x_j - x_i >= 0
                    data[pos_i] = -1.0;
                    data[pos_j] = 1.0;
                } else {
                    // x_i - x_j >= 0
                    data[pos_i] = 1.0;
                    data[pos_j] = -1.0;
                }
            }

            // Optional: reset marker (not strictly necessary since each row should be hit exactly twice).
            // first_col[r] = usize::MAX;
        }
    }

    debug_assert!(
        first_col.iter().all(|&c| c != usize::MAX),
        "some rows did not have exactly two nonzeros"
    );
}

#[cfg(test)]
mod tests {
    use super::update_constraint_matrix;
    use crate::partial_order::algorithm::uncensored::build_order_constraints;
    use crate::structures::{Decreasing, Increasing};
    use osqp::CscMatrix;

    // Find the position in data/indices for (col, row).
    fn pos(m: &CscMatrix, col: usize, row: usize) -> usize {
        let start = m.indptr[col];
        let end = m.indptr[col + 1];
        for k in start..end {
            if m.indices[k] == row {
                return k;
            }
        }
        panic!("Row {row} not found in column {col}");
    }

    // Build the expected data vector after an update for a given survival vector.
    //
    // The rule under test:
    // - If S[i_row] == 0, row r = (i, j) uses the "simple" builder signs (base_simple).
    // - If S[i_row] > 0, row uses the "hazard" builder signs (base_hazard), scaled as:
    //     value at i becomes (base_hazard_sign at i)/S_i,
    //     value at j becomes (base_hazard_sign at j)/S_j.
    //
    // This definition works whether the current matrix before update is simple or hazard
    // (the function must map from either state to this canonical target).
    fn expected_data_from_bases(
        base_simple: &CscMatrix,
        base_hazard: &CscMatrix,
        constraints: &[(usize, usize)],
        survival: &[f64],
    ) -> Vec<f64> {
        let mut expected = base_simple.data.to_vec(); // start with simple; we’ll overwrite where S_i > 0.

        for (r, &(i, j)) in constraints.iter().enumerate() {
            let ki = pos(base_simple, i, r);
            let kj = pos(base_simple, j, r);
            if survival[i] == 0.0 {
                // Keep simple signs (already in expected from base_simple).
                expected[ki] = base_simple.data[ki];
                expected[kj] = base_simple.data[kj];
            } else {
                // Use hazard signs scaled by reciprocals.
                let ki_h = pos(base_hazard, i, r);
                let kj_h = pos(base_hazard, j, r);
                expected[ki] = base_hazard.data[ki_h] / survival[i];
                expected[kj] = base_hazard.data[kj_h] / survival[j];
            }
        }
        expected
    }

    fn assert_same_structure(a: &CscMatrix, b: &CscMatrix) {
        assert_eq!(a.nrows, b.nrows, "nrows changed");
        assert_eq!(a.ncols, b.ncols, "ncols changed");
        assert_eq!(a.indptr.as_ref(), b.indptr.as_ref(), "indptr changed");
        assert_eq!(a.indices.as_ref(), b.indices.as_ref(), "indices changed");
    }

    fn assert_data_eq(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len(), "data length mismatch");
        for (k, (a, e)) in actual
            .iter()
            .copied()
            .zip(expected.iter().copied())
            .enumerate()
        {
            // We pick survival values so these are exact binary floats, but allow a tiny epsilon anyway.
            let eps = 1e-12;
            assert!((a - e).abs() <= eps, "data[{k}] = {a} != {e}");
        }
    }

    #[test]
    fn structure_is_preserved() {
        // Simple 1-edge case.
        let constraints = vec![(0usize, 1usize)];
        let base_simple = build_order_constraints::<Increasing, false>(&constraints, 2);
        let mut mat = base_simple.clone();
        let survival = vec![0.5, 1.0];

        update_constraint_matrix::<Increasing>(&mut mat, &survival);

        assert_same_structure(&mat, &base_simple);
        // Two nonzeros total; and two per row in this case.
        assert_eq!(mat.indptr[0], 0);
        assert_eq!(mat.indptr[1] - mat.indptr[0], 1);
        assert_eq!(mat.indptr[2] - mat.indptr[1], 1);
    }

    #[test]
    fn increasing_basic_no_zeros_and_idempotent() {
        // Graph: 0 -> 1 -> 2
        let constraints = vec![(0, 1), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Increasing, true>(&constraints, n);

        // survival satisfies S0 ≤ S1 ≤ S2
        let survival = vec![0.5, 0.75, 1.0];

        let mut mat = base_simple.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &survival);

        // Expected: hazard signs scaled by reciprocals.
        let expected =
            expected_data_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_same_structure(&mat, &base_simple);
        assert_data_eq(mat.data.as_ref(), &expected);

        // Idempotence: applying again with same S yields identical data.
        let before = mat.data.to_owned().into_owned();
        update_constraint_matrix::<Increasing>(&mut mat, &survival);
        assert_data_eq(mat.data.as_ref(), &before.as_slice());
    }

    #[test]
    fn increasing_zero_at_i_reverts_row_to_simple() {
        // Graph: 0 -> 1 -> 2 (only row 0 uses i=0)
        let constraints = vec![(0, 1), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Increasing, true>(&constraints, n);

        // Make S0 = 0 (and maintain monotone: 0 ≤ 0.5 ≤ 1.0).
        let survival = vec![0.0, 0.5, 1.0];

        let mut mat = base_simple.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &survival);

        let expected =
            expected_data_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_same_structure(&mat, &base_simple);
        assert_data_eq(mat.data.as_ref(), &expected);
    }

    #[test]
    fn increasing_multiple_edges_shared_columns() {
        // Graph: edges (0,1), (0,2), (1,2)
        let constraints = vec![(0, 1), (0, 2), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Increasing, true>(&constraints, n);

        // survival: 0.25 ≤ 0.5 ≤ 1.0
        let survival = vec![0.25, 0.5, 1.0];

        let mut mat = base_simple.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &survival);

        let expected =
            expected_data_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_same_structure(&mat, &base_simple);
        assert_data_eq(mat.data.as_ref(), &expected);
    }

    #[test]
    fn increasing_iterative_updates_with_changing_survival() {
        let constraints = vec![(0, 1), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Increasing, true>(&constraints, n);

        let mut mat = base_simple.clone();

        // First survival (strictly positive).
        let s1 = vec![0.5, 0.75, 1.0];
        update_constraint_matrix::<Increasing>(&mut mat, &s1);
        let expected1 = expected_data_from_bases(&base_simple, &base_hazard, &constraints, &s1);
        assert_data_eq(mat.data.as_ref(), &expected1);

        // Change survival: set S1 to 0 (row 1 should revert), others change scaling.
        let s2 = vec![0.0, 0.5, 1.0];
        update_constraint_matrix::<Increasing>(&mut mat, &s2);
        let expected2 = expected_data_from_bases(&base_simple, &base_hazard, &constraints, &s2);
        assert_data_eq(mat.data.as_ref(), &expected2);

        // Back to positive S1; should return to hazard scaling.
        let s3 = vec![0.5, 0.5, 1.0];
        update_constraint_matrix::<Increasing>(&mut mat, &s3);
        let expected3 = expected_data_from_bases(&base_simple, &base_hazard, &constraints, &s3);
        assert_data_eq(mat.data.as_ref(), &expected3);
    }

    #[test]
    fn starting_from_hazard_oriented_builder() {
        // Same constraints, but the matrix was created with HRO=true initially.
        let constraints = vec![(0, 1), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Increasing, true>(&constraints, n);

        // survival strictly positive
        let survival = vec![0.25, 0.5, 1.0];

        // Start from hazard-oriented matrix.
        let mut mat = base_hazard.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &survival);

        // Expected does not depend on the current state; it must be "hazard signs scaled", or
        // "simple signs" for zero S_i (there are none here).
        let expected =
            expected_data_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_same_structure(&mat, &base_hazard);
        assert_data_eq(mat.data.as_ref(), &expected);
    }

    #[test]
    fn decreasing_direction_all_cases_mirror() {
        // Start from a matrix built with HRO = true (this is the case we care about).
        let constraints = vec![(0usize, 1usize), (1usize, 2usize)];
        let n = 3usize;

        // Build the base (hazard-oriented) matrix for decreasing direction.
        let base_hazard = build_order_constraints::<Decreasing, true>(&constraints, n);
        let mut mat = base_hazard.clone();

        // Helper: position in CSC storage for (col, row), using the (unchanging) base structure.
        let pos = |col: usize, row: usize| -> usize {
            let start = base_hazard.indptr[col];
            let end = base_hazard.indptr[col + 1];
            for k in start..end {
                if base_hazard.indices[k] == row {
                    return k;
                }
            }
            panic!("Row {row} not found in column {col}");
        };

        // Reproduce the implementation’s orientation and coefficient rule for Decreasing:
        // - i = endpoint with larger (or equal) survival (tie-break: larger column index).
        // - If both survivals > 0: set -1/S_i at i and +1/S_j at j.
        // - Else: fallback simple decreasing: +1 at i and -1 at j.
        let mat_nrows = mat.nrows;
        let expected_for = |survival: &[f64]| -> Vec<f64> {
            let mrows = mat_nrows;
            let nnz = base_hazard.data.len();
            let mut out = vec![0.0f64; nnz];

            for r in 0..mrows {
                let (c0, c1) = constraints[r];
                let s0 = survival[c0];
                let s1 = survival[c1];

                // orient_pair for decreasing (match the implementation)
                let (i_col, s_i, j_col, s_j) = if s0 > s1 || (s0 == s1 && c0 >= c1) {
                    (c0, s0, c1, s1)
                } else {
                    (c1, s1, c0, s0)
                };

                let pi = pos(i_col, r);
                let pj = pos(j_col, r);

                if s_i > 0.0 && s_j > 0.0 {
                    out[pi] = -1.0 / s_i;
                    out[pj] = 1.0 / s_j;
                } else {
                    out[pi] = 1.0;
                    out[pj] = -1.0;
                }
            }
            out
        };

        // Case 1: strictly positive survivals (S0 ≥ S1 ≥ S2).
        let s1 = vec![1.0, 0.5, 0.25];
        update_constraint_matrix::<Decreasing>(&mut mat, &s1);

        // Structure must be preserved.
        assert_eq!(mat.nrows, base_hazard.nrows);
        assert_eq!(mat.ncols, base_hazard.ncols);
        assert_eq!(mat.indptr.as_ref(), base_hazard.indptr.as_ref());
        assert_eq!(mat.indices.as_ref(), base_hazard.indices.as_ref());

        // Data must match the hazard-scaled coefficients.
        let expected1 = expected_for(&s1);
        for (k, (a, e)) in mat
            .data
            .iter()
            .copied()
            .zip(expected1.iter().copied())
            .enumerate()
        {
            assert!((a - e).abs() <= 1e-12, "after s1: data[{k}] = {a} != {e}");
        }

        // Case 2: introduce zeros to trigger fallback on rows touching zero-survival endpoints.
        // Here S1 = 0 and S2 = 0; orientation for row (1,2) will pick i=2, j=1 due to tie-break.
        let s2 = vec![1.0, 0.0, 0.0];
        super::update_constraint_matrix::<Decreasing>(&mut mat, &s2);

        let expected2 = expected_for(&s2);
        for (k, (a, e)) in mat
            .data
            .iter()
            .copied()
            .zip(expected2.iter().copied())
            .enumerate()
        {
            assert!((a - e).abs() <= 1e-12, "after s2: data[{k}] = {a} != {e}");
        }
    }
}
