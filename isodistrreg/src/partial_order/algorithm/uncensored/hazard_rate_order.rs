use crate::structures::Direction;
use osqp::CscMatrix;
use std::borrow::Cow;

#[allow(clippy::doc_overindented_list_items)]
/// Update A in-place to encode HRO constraints per edge-row.
///
/// The QP runs in energy coordinates `u = sqrt_weight[c] * x_c` (see `algorithm`), where
/// the x-variables are the survival values at the current threshold and `S` holds the
/// survivals at the previous threshold, so `x/S` is the per-step survival ratio. The
/// hazard rate order requires that ratio to be nondecreasing along the covariate order
/// (Y_i <=hr Y_j iff S_j(t)/S_i(t) is nondecreasing in t) — the same orientation as the
/// plain-order fallback, just ratio-weighted.
///
/// `constraints[r] = (i, j)` is the transitive-reduction edge encoded by row r, with i
/// the covariate that is below j in the covariate order. The orientation of each row is
/// fully determined by that statically known edge direction — it must NOT be inferred
/// from the previous fitted survivals: pooled blocks make those survivals tied up to
/// solver round-off, which would let noise decide the constraint direction.
///
/// For each row r encoding the edge (i, j), in x-space:
/// - If D::IS_INCREASING:
///     HRO:  m * ((x_j / S_j) - (x_i / S_i)) >= 0, if S_i > 0 and S_j > 0
///     else: m' *             (x_j - x_i)    >= 0
/// - Else (decreasing):
///     HRO:  m * ((x_i / S_i) - (x_j / S_j)) >= 0, if S_i > 0 and S_j > 0
///     else: m' *             (x_i - x_j)    >= 0
///
/// In energy coordinates the column-`c` coefficient additionally divides by
/// `sqrt_weight[c]`, so with `d_c = sqrt_weight[c] * S_c` the stored coefficients are
/// ±m/d_i and ∓m/d_j with the row scaling `m = min(d_i, d_j)` (and, for the plain-order
/// fallback, ±m'/sqrt_weight with `m' = min(sqrt_weight[i], sqrt_weight[j])`). The
/// positive row scaling keeps the constraint equivalent while bounding both coefficients
/// by 1 in magnitude, so the row value stays in [-1, 1] for x in [0, 1]^n and the QP's
/// fixed upper bound of 1.0 per row (an OSQP duality-gap workaround, see `algorithm`)
/// remains redundant. Unscaled coefficients can produce row values above 1 that the
/// bound would clip.
///
/// Invariants expected:
/// - existing_matrix has exactly two nonzeros per row, in the columns given by
///   `constraints` (built by `build_order_constraints` or previously modified only by
///   this function or `scale_constraints_to_energy`).
/// - constraints.len() == existing_matrix.nrows.
/// - survival.len() == sqrt_weight.len() == existing_matrix.ncols.
/// - sqrt_weight is strictly positive.
pub fn update_constraint_matrix<D: Direction>(
    existing_matrix: &mut CscMatrix,
    constraints: &[(usize, usize)],
    survival: &[f64],
    sqrt_weight: &[f64],
) {
    let m = existing_matrix.nrows;
    let n = existing_matrix.ncols;

    assert_eq!(
        constraints.len(),
        m,
        "constraints length ({}) must equal number of rows ({})",
        constraints.len(),
        m
    );
    assert_eq!(
        survival.len(),
        n,
        "survival length ({}) must equal number of columns ({})",
        survival.len(),
        n
    );
    assert_eq!(
        sqrt_weight.len(),
        n,
        "sqrt_weight length ({}) must equal number of columns ({})",
        sqrt_weight.len(),
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

    // Locate the storage position of each row's (i, j) entries. usize::MAX = not found.
    let mut positions = vec![(usize::MAX, usize::MAX); m];
    for col in 0..n {
        #[allow(clippy::needless_range_loop)]
        for pos in indptr[col]..indptr[col + 1] {
            let r = indices[pos];
            debug_assert!(r < m, "row index out of range");
            let (i, j) = constraints[r];
            if col == i {
                positions[r].0 = pos;
            } else {
                debug_assert_eq!(col, j, "row {r} has an entry outside its edge columns");
                positions[r].1 = pos;
            }
        }
    }

    for (r, &(i, j)) in constraints.iter().enumerate() {
        let (pos_i, pos_j) = positions[r];
        debug_assert!(
            pos_i != usize::MAX && pos_j != usize::MAX,
            "row {r} does not have exactly two nonzeros"
        );

        let s_i = survival[i];
        let s_j = survival[j];

        if s_i > 0.0 && s_j > 0.0 {
            // `d_c` is the u-extent of column c: the ratio constraint divides x_c by S_c
            // and energy coordinates divide by sqrt_weight[c]. The equivalent positive row
            // rescaling by min(d_i, d_j) keeps |coefficient| <= 1.
            let d_i = sqrt_weight[i] * s_i;
            let d_j = sqrt_weight[j] * s_j;
            let scale = d_i.min(d_j);
            if D::IS_INCREASING {
                // (x_j/S_j) - (x_i/S_i) >= 0
                data[pos_i] = -scale / d_i;
                data[pos_j] = scale / d_j;
            } else {
                // (x_i/S_i) - (x_j/S_j) >= 0
                data[pos_i] = scale / d_i;
                data[pos_j] = -scale / d_j;
            }
        } else {
            // Fallback to simple order (in energy coordinates, with the same row scaling
            // rule as `scale_constraints_to_energy`).
            let scale = sqrt_weight[i].min(sqrt_weight[j]);
            if D::IS_INCREASING {
                // x_j - x_i >= 0
                data[pos_i] = -scale / sqrt_weight[i];
                data[pos_j] = scale / sqrt_weight[j];
            } else {
                // x_i - x_j >= 0
                data[pos_i] = scale / sqrt_weight[i];
                data[pos_j] = -scale / sqrt_weight[j];
            }
        }
    }
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

    // Build the expected data vector after an update for a given survival vector, with
    // unit sqrt-weights (energy coordinates with sqrt_weight = 1 reduce to x-space; see
    // `energy_coordinates_scale_columns` for the non-unit rule).
    //
    // The rule under test (orientation is static, from the row's edge (i, j)):
    // - If S_i > 0 and S_j > 0, row r = (i, j) uses the "hazard" builder signs
    //   (base_hazard), scaled with the normalization m = min(S_i, S_j):
    //     value at i becomes (base_hazard_sign at i) * m / S_i,
    //     value at j becomes (base_hazard_sign at j) * m / S_j.
    // - Otherwise the row falls back to the "simple" builder signs (base_simple).
    //
    // This definition works whether the current matrix before update is simple or hazard
    // (the function must map from either state to this canonical target).
    fn expected_data_from_bases(
        base_simple: &CscMatrix,
        base_hazard: &CscMatrix,
        constraints: &[(usize, usize)],
        survival: &[f64],
    ) -> Vec<f64> {
        let mut expected = base_simple.data.to_vec(); // start with simple; we’ll overwrite where S > 0.

        for (r, &(i, j)) in constraints.iter().enumerate() {
            let ki = pos(base_simple, i, r);
            let kj = pos(base_simple, j, r);
            if survival[i] > 0.0 && survival[j] > 0.0 {
                // Use hazard signs scaled by normalized reciprocals.
                let scale = survival[i].min(survival[j]);
                let ki_h = pos(base_hazard, i, r);
                let kj_h = pos(base_hazard, j, r);
                expected[ki] = base_hazard.data[ki_h] * scale / survival[i];
                expected[kj] = base_hazard.data[kj_h] * scale / survival[j];
            } else {
                // Keep simple signs (already in expected from base_simple).
                expected[ki] = base_simple.data[ki];
                expected[kj] = base_simple.data[kj];
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

        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &[1.0; 2]);

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
        let base_hazard = build_order_constraints::<Decreasing, true>(&constraints, n);

        // survival satisfies S0 ≤ S1 ≤ S2
        let survival = vec![0.5, 0.75, 1.0];

        let mut mat = base_simple.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &[1.0; 3]);

        // Expected: hazard signs scaled by reciprocals.
        let expected =
            expected_data_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_same_structure(&mat, &base_simple);
        assert_data_eq(mat.data.as_ref(), &expected);

        // Idempotence: applying again with same S yields identical data.
        let before = mat.data.to_owned().into_owned();
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &[1.0; 3]);
        assert_data_eq(mat.data.as_ref(), &before.as_slice());
    }

    #[test]
    fn increasing_zero_at_i_reverts_row_to_simple() {
        // Graph: 0 -> 1 -> 2 (only row 0 uses i=0)
        let constraints = vec![(0, 1), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Decreasing, true>(&constraints, n);

        // Make S0 = 0 (and maintain monotone: 0 ≤ 0.5 ≤ 1.0).
        let survival = vec![0.0, 0.5, 1.0];

        let mut mat = base_simple.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &[1.0; 3]);

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
        let base_hazard = build_order_constraints::<Decreasing, true>(&constraints, n);

        // survival: 0.25 ≤ 0.5 ≤ 1.0
        let survival = vec![0.25, 0.5, 1.0];

        let mut mat = base_simple.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &[1.0; 3]);

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
        let base_hazard = build_order_constraints::<Decreasing, true>(&constraints, n);

        let mut mat = base_simple.clone();

        // First survival (strictly positive).
        let s1 = vec![0.5, 0.75, 1.0];
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &s1, &[1.0; 3]);
        let expected1 = expected_data_from_bases(&base_simple, &base_hazard, &constraints, &s1);
        assert_data_eq(mat.data.as_ref(), &expected1);

        // Change survival: set S1 to 0 (row 1 should revert), others change scaling.
        let s2 = vec![0.0, 0.5, 1.0];
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &s2, &[1.0; 3]);
        let expected2 = expected_data_from_bases(&base_simple, &base_hazard, &constraints, &s2);
        assert_data_eq(mat.data.as_ref(), &expected2);

        // Back to positive S1; should return to hazard scaling.
        let s3 = vec![0.5, 0.5, 1.0];
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &s3, &[1.0; 3]);
        let expected3 = expected_data_from_bases(&base_simple, &base_hazard, &constraints, &s3);
        assert_data_eq(mat.data.as_ref(), &expected3);
    }

    /// Non-unit sqrt-weights (energy coordinates): column c divides by
    /// sqrt_weight[c] * S_c, and the row is rescaled so the largest |coefficient| is
    /// exactly 1. All values below are exact binary floats.
    #[test]
    fn energy_coordinates_scale_columns() {
        let constraints = vec![(0usize, 1usize)];
        let base_simple = build_order_constraints::<Increasing, false>(&constraints, 2);
        let mut mat = base_simple.clone();
        let sqrt_weight = [1.0, 0.125];

        // Both survivals positive: d_0 = 1.0 * 0.5 = 0.5, d_1 = 0.125 * 0.25 = 0.03125,
        // row scale min(d) = 0.03125 -> -0.0625 at column 0 and +1.0 at column 1.
        let survival = [0.5, 0.25];
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &sqrt_weight);
        assert_same_structure(&mat, &base_simple);
        let k0 = pos(&base_simple, 0, 0);
        let k1 = pos(&base_simple, 1, 0);
        assert_eq!(mat.data[k0], -0.0625);
        assert_eq!(mat.data[k1], 1.0);

        // Zero survival falls back to the plain order in energy coordinates:
        // row scale min(sqrt_weight) = 0.125 -> -0.125 at column 0 and +1.0 at column 1.
        let survival = [0.5, 0.0];
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &sqrt_weight);
        assert_eq!(mat.data[k0], -0.125);
        assert_eq!(mat.data[k1], 1.0);
    }

    #[test]
    fn starting_from_hazard_oriented_builder() {
        // Same constraints, but the matrix was created with HRO=true initially.
        let constraints = vec![(0, 1), (1, 2)];
        let n = 3;

        let base_simple = build_order_constraints::<Increasing, false>(&constraints, n);
        let base_hazard = build_order_constraints::<Decreasing, true>(&constraints, n);

        // survival strictly positive
        let survival = vec![0.25, 0.5, 1.0];

        // Start from hazard-oriented matrix.
        let mut mat = base_hazard.clone();
        update_constraint_matrix::<Increasing>(&mut mat, &constraints, &survival, &[1.0; 3]);

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

        // Build the base (hazard-oriented) matrix for a decreasing fit — the pipeline
        // builds with D::REVERSE (algorithm/uncensored/mod.rs), so Increasing here.
        let base_hazard = build_order_constraints::<Increasing, true>(&constraints, n);
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

        // Reproduce the coefficient rule for Decreasing, with i the edge's statically
        // known lower endpoint (constraints[r].0):
        // - If both survivals > 0: with m = min(S_i, S_j), set +m/S_i at i and -m/S_j
        //   at j, i.e. m * ((x_i/S_i) - (x_j/S_j)) >= 0 — the survival ratio
        //   nonincreasing along the covariate order, the hazard rate order for a
        //   decreasing fit, rescaled so the row's coefficients stay within ±1.
        // - Else: fallback simple decreasing: +1 at i and -1 at j.
        let mat_nrows = mat.nrows;
        let expected_for = |survival: &[f64]| -> Vec<f64> {
            let mrows = mat_nrows;
            let nnz = base_hazard.data.len();
            let mut out = vec![0.0f64; nnz];

            for r in 0..mrows {
                let (i_col, j_col) = constraints[r];
                let s_i = survival[i_col];
                let s_j = survival[j_col];

                let pi = pos(i_col, r);
                let pj = pos(j_col, r);

                if s_i > 0.0 && s_j > 0.0 {
                    let scale = s_i.min(s_j);
                    out[pi] = scale / s_i;
                    out[pj] = -scale / s_j;
                } else {
                    out[pi] = 1.0;
                    out[pj] = -1.0;
                }
            }
            out
        };

        // Case 1: strictly positive survivals (S0 ≥ S1 ≥ S2).
        let s1 = vec![1.0, 0.5, 0.25];
        update_constraint_matrix::<Decreasing>(&mut mat, &constraints, &s1, &[1.0; 3]);

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
        let s2 = vec![1.0, 0.0, 0.0];
        super::update_constraint_matrix::<Decreasing>(&mut mat, &constraints, &s2, &[1.0; 3]);

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
