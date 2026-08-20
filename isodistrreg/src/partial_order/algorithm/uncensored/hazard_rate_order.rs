use crate::structures::Direction;
/// Rewrite the constraint coefficients for the current threshold's ratio constraints.
///
/// See the module documentation for the orientation rules. `coefficients[r]` is the pair
/// `(coefficient at i, coefficient at j)` for the edge `constraints[r] = (i, j)`; the
/// write is absolute, not multiplicative, so this maps coefficients in any previous state
/// to the canonical target and is idempotent.
pub fn update_constraint_coefficients<D: Direction>(
    coefficients: &mut [(f64, f64)],
    constraints: &[(usize, usize)],
    survival: &[f64],
    sqrt_weight: &[f64],
) {
    debug_assert_eq!(coefficients.len(), constraints.len());
    debug_assert_eq!(survival.len(), sqrt_weight.len());

    for (coefficient, &(i, j)) in coefficients.iter_mut().zip(constraints) {
        let s_i = survival[i];
        let s_j = survival[j];

        *coefficient = if s_i > 0.0 && s_j > 0.0 {
            // `d_c` is the u-extent of column c: the ratio constraint divides x_c by S_c
            // and energy coordinates divide by sqrt_weight[c]. The equivalent positive row
            // rescaling by min(d_i, d_j) keeps |coefficient| <= 1.
            let d_i = sqrt_weight[i] * s_i;
            let d_j = sqrt_weight[j] * s_j;
            let scale = d_i.min(d_j);
            if D::IS_INCREASING {
                // (x_j/S_j) - (x_i/S_i) >= 0
                (-scale / d_i, scale / d_j)
            } else {
                // (x_i/S_i) - (x_j/S_j) >= 0
                (scale / d_i, -scale / d_j)
            }
        } else {
            // Fallback to simple order (in energy coordinates, with the same row scaling
            // rule as `scale_coefficients_to_energy`).
            let scale = sqrt_weight[i].min(sqrt_weight[j]);
            if D::IS_INCREASING {
                // x_j - x_i >= 0
                (-scale / sqrt_weight[i], scale / sqrt_weight[j])
            } else {
                // x_i - x_j >= 0
                (scale / sqrt_weight[i], -scale / sqrt_weight[j])
            }
        };
    }
}

#[cfg(test)]
mod tests {
    use super::update_constraint_coefficients;
    use crate::partial_order::algorithm::uncensored::build_order_coefficients;
    use crate::structures::{Decreasing, Increasing};

    /// The canonical coefficients an update must produce, whatever state it starts from.
    ///
    /// Orientation is static, taken from the row's edge `(i, j)`:
    /// - if `S_i > 0` and `S_j > 0`, the row uses the hazard-oriented builder's signs
    ///   scaled by `m = min(S_i, S_j)`: `sign_i * m / S_i` at `i`, `sign_j * m / S_j` at `j`;
    /// - otherwise it falls back to the plain-order builder's signs.
    ///
    /// Unit sqrt-weights, so energy coordinates reduce to x-space; the non-unit rule is
    /// checked separately by `energy_coordinates_scale_columns`.
    fn expected_from_bases(
        base_simple: &[(f64, f64)],
        base_hazard: &[(f64, f64)],
        constraints: &[(usize, usize)],
        survival: &[f64],
    ) -> Vec<(f64, f64)> {
        constraints
            .iter()
            .enumerate()
            .map(|(r, &(i, j))| {
                if survival[i] > 0.0 && survival[j] > 0.0 {
                    let scale = survival[i].min(survival[j]);
                    (
                        base_hazard[r].0 * scale / survival[i],
                        base_hazard[r].1 * scale / survival[j],
                    )
                } else {
                    base_simple[r]
                }
            })
            .collect()
    }

    /// Survival values are chosen so the expected coefficients are exact binary floats,
    /// but a tiny epsilon is allowed anyway.
    fn assert_coefficients_eq(actual: &[(f64, f64)], expected: &[(f64, f64)]) {
        assert_eq!(actual.len(), expected.len(), "row count mismatch");
        for (r, (a, e)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (a.0 - e.0).abs() <= 1e-12 && (a.1 - e.1).abs() <= 1e-12,
                "row {r}: {a:?} != {e:?}"
            );
        }
    }

    /// The row-to-edge pairing survives an update: one coefficient per endpoint, opposite
    /// signs, in the order the edge lists them. With a CSC matrix this needed an explicit
    /// structure check; here the edge index *is* the storage index, so what remains to
    /// verify is that the update writes both endpoints and preserves the orientation.
    #[test]
    fn row_pairing_is_preserved() {
        let constraints = vec![(0usize, 1usize)];
        let mut coefficients = build_order_coefficients::<Increasing, false>(&constraints);
        let survival = vec![0.5, 1.0];

        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &[1.0; 2],
        );

        assert_eq!(coefficients.len(), constraints.len());
        let (at_i, at_j) = coefficients[0];
        assert!(
            at_i < 0.0 && at_j > 0.0,
            "orientation lost: {:?}",
            coefficients[0]
        );
    }

    #[test]
    fn increasing_basic_no_zeros_and_idempotent() {
        // Graph: 0 -> 1 -> 2
        let constraints = vec![(0, 1), (1, 2)];

        let base_simple = build_order_coefficients::<Increasing, false>(&constraints);
        let base_hazard = build_order_coefficients::<Decreasing, true>(&constraints);

        // survival satisfies S0 <= S1 <= S2
        let survival = vec![0.5, 0.75, 1.0];

        let mut coefficients = base_simple.clone();
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &[1.0; 3],
        );

        let expected = expected_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_coefficients_eq(&coefficients, &expected);

        // Idempotence: applying again with the same S yields identical coefficients.
        let before = coefficients.clone();
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &[1.0; 3],
        );
        assert_coefficients_eq(&coefficients, &before);
    }

    #[test]
    fn increasing_zero_at_i_reverts_row_to_simple() {
        // Graph: 0 -> 1 -> 2 (only row 0 uses i = 0)
        let constraints = vec![(0, 1), (1, 2)];

        let base_simple = build_order_coefficients::<Increasing, false>(&constraints);
        let base_hazard = build_order_coefficients::<Decreasing, true>(&constraints);

        // Make S0 = 0 (and maintain monotone: 0 <= 0.5 <= 1.0).
        let survival = vec![0.0, 0.5, 1.0];

        let mut coefficients = base_simple.clone();
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &[1.0; 3],
        );

        let expected = expected_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_coefficients_eq(&coefficients, &expected);
    }

    #[test]
    fn increasing_multiple_edges_shared_columns() {
        // Graph: edges (0,1), (0,2), (1,2)
        let constraints = vec![(0, 1), (0, 2), (1, 2)];

        let base_simple = build_order_coefficients::<Increasing, false>(&constraints);
        let base_hazard = build_order_coefficients::<Decreasing, true>(&constraints);

        // survival: 0.25 <= 0.5 <= 1.0
        let survival = vec![0.25, 0.5, 1.0];

        let mut coefficients = base_simple.clone();
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &[1.0; 3],
        );

        let expected = expected_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_coefficients_eq(&coefficients, &expected);
    }

    #[test]
    fn increasing_iterative_updates_with_changing_survival() {
        let constraints = vec![(0, 1), (1, 2)];

        let base_simple = build_order_coefficients::<Increasing, false>(&constraints);
        let base_hazard = build_order_coefficients::<Decreasing, true>(&constraints);

        let mut coefficients = base_simple.clone();

        // First survival (strictly positive).
        let s1 = vec![0.5, 0.75, 1.0];
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &s1,
            &[1.0; 3],
        );
        assert_coefficients_eq(
            &coefficients,
            &expected_from_bases(&base_simple, &base_hazard, &constraints, &s1),
        );

        // Change survival: set S0 to 0 (row 0 should revert), others change scaling.
        let s2 = vec![0.0, 0.5, 1.0];
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &s2,
            &[1.0; 3],
        );
        assert_coefficients_eq(
            &coefficients,
            &expected_from_bases(&base_simple, &base_hazard, &constraints, &s2),
        );

        // Back to positive S0; should return to hazard scaling.
        let s3 = vec![0.5, 0.5, 1.0];
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &s3,
            &[1.0; 3],
        );
        assert_coefficients_eq(
            &coefficients,
            &expected_from_bases(&base_simple, &base_hazard, &constraints, &s3),
        );
    }

    /// Non-unit sqrt-weights (energy coordinates): column c divides by
    /// `sqrt_weight[c] * S_c`, and the row is rescaled so the largest |coefficient| is
    /// exactly 1. All values below are exact binary floats.
    #[test]
    fn energy_coordinates_scale_columns() {
        let constraints = vec![(0usize, 1usize)];
        let mut coefficients = build_order_coefficients::<Increasing, false>(&constraints);
        let sqrt_weight = [1.0, 0.125];

        // Both survivals positive: d_0 = 1.0 * 0.5 = 0.5, d_1 = 0.125 * 0.25 = 0.03125,
        // row scale min(d) = 0.03125 -> -0.0625 at column 0 and +1.0 at column 1.
        let survival = [0.5, 0.25];
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &sqrt_weight,
        );
        assert_eq!(coefficients[0], (-0.0625, 1.0));

        // Zero survival falls back to the plain order in energy coordinates:
        // row scale min(sqrt_weight) = 0.125 -> -0.125 at column 0 and +1.0 at column 1.
        let survival = [0.5, 0.0];
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &sqrt_weight,
        );
        assert_eq!(coefficients[0], (-0.125, 1.0));
    }

    #[test]
    fn starting_from_hazard_oriented_builder() {
        // Same constraints, but the coefficients were created with HRO = true initially.
        let constraints = vec![(0, 1), (1, 2)];

        let base_simple = build_order_coefficients::<Increasing, false>(&constraints);
        let base_hazard = build_order_coefficients::<Decreasing, true>(&constraints);

        // survival strictly positive
        let survival = vec![0.25, 0.5, 1.0];

        // Start from the hazard-oriented coefficients.
        let mut coefficients = base_hazard.clone();
        update_constraint_coefficients::<Increasing>(
            &mut coefficients,
            &constraints,
            &survival,
            &[1.0; 3],
        );

        // The target does not depend on the current state.
        let expected = expected_from_bases(&base_simple, &base_hazard, &constraints, &survival);
        assert_coefficients_eq(&coefficients, &expected);
    }

    #[test]
    fn decreasing_direction_all_cases_mirror() {
        // Start from coefficients built with HRO = true (this is the case we care about).
        let constraints = vec![(0usize, 1usize), (1usize, 2usize)];

        // The pipeline builds with D::REVERSE (algorithm/uncensored/mod.rs), so a
        // decreasing fit builds its base with Increasing here.
        let base_hazard = build_order_coefficients::<Increasing, true>(&constraints);
        let mut coefficients = base_hazard.clone();

        // Reproduce the coefficient rule for Decreasing, with i the edge's statically
        // known lower endpoint (constraints[r].0):
        // - If both survivals > 0: with m = min(S_i, S_j), set +m/S_i at i and -m/S_j
        //   at j, i.e. m * ((x_i/S_i) - (x_j/S_j)) >= 0 -- the survival ratio
        //   nonincreasing along the covariate order, the hazard rate order for a
        //   decreasing fit, rescaled so the row's coefficients stay within +-1.
        // - Else: fallback simple decreasing: +1 at i and -1 at j.
        let expected_for = |survival: &[f64]| -> Vec<(f64, f64)> {
            constraints
                .iter()
                .map(|&(i, j)| {
                    let (s_i, s_j) = (survival[i], survival[j]);
                    if s_i > 0.0 && s_j > 0.0 {
                        let scale = s_i.min(s_j);
                        (scale / s_i, -scale / s_j)
                    } else {
                        (1.0, -1.0)
                    }
                })
                .collect()
        };

        // Case 1: strictly positive survivals (S0 >= S1 >= S2).
        let s1 = vec![1.0, 0.5, 0.25];
        update_constraint_coefficients::<Decreasing>(
            &mut coefficients,
            &constraints,
            &s1,
            &[1.0; 3],
        );
        assert_coefficients_eq(&coefficients, &expected_for(&s1));

        // Case 2: introduce zeros to trigger fallback on rows touching zero-survival
        // endpoints.
        let s2 = vec![1.0, 0.0, 0.0];
        update_constraint_coefficients::<Decreasing>(
            &mut coefficients,
            &constraints,
            &s2,
            &[1.0; 3],
        );
        assert_coefficients_eq(&coefficients, &expected_for(&s2));
    }
}
