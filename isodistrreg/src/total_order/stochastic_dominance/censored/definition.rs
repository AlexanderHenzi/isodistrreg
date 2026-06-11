//! Total-order reference implementation: the chain specialization of the general
//! partial-order estimator
//!
//! ```text
//! F̂_{ξ_i}(y) = min_{U ∈ 𝓤: i ∈ U}  max_{L ∈ 𝓛: i ∈ L}  RKM_{L ∩ U}(y)
//! RKM_S(y)   = clamp( KM_S(y),
//!                     max_{(L,U) ∈ P_S} min(RKM_{S∩L}, RKM_{S∩U}),
//!                     min_{(L,U) ∈ P_S} max(RKM_{S∩L}, RKM_{S∩U}) )
//! P_S        = { (S∩L, S∩U) : L ∈ 𝓛, U ∈ 𝓤, L ∪ U ⊇ S, L ∩ U = ∅ }
//! ```
//!
//! (see `partial_order::algorithm::definition` for the literal poset implementation).
//! On a totally ordered covariate grid the lower sets 𝓛 are the prefixes, the upper
//! sets 𝓤 the suffixes, every L ∩ U is an interval `[r, s]`, and the nondegenerate
//! pairs of P_{[r,s]} are exactly the split points `k ∈ r..s` — the degenerate pairs
//! (S, ∅)/(∅, S) are excluded since RKM_∅ is undefined, here implicitly, in the poset
//! implementations explicitly.
//!
//! This module computes in SURVIVAL space, where `1 − x` swaps min ↔ max:
//! - the clip bounds `max_k min(..) / min_k max(..)` below are the `1 − x` image of
//!   the spec's clamp bounds;
//! - clipping in span-ascending order makes every clip input an already-clipped value,
//!   i.e. the recursion on RKM (not raw KM) values, mirrored by popcount-ascending
//!   memoization in the poset spec;
//! - the outer reduction (Increasing: inner `min` over interval ends `s ≥ i` = the
//!   CDF-space `max_L` over prefixes, outer `max` over interval starts `r ≤ i` = the
//!   CDF-space `min_U` over suffixes) matches the spec's nesting order exactly — the
//!   correspondence is syntactic, no minimax/saddle interchange is involved.
//!
//! The correspondence is enforced by the differential suites in `super::test` (the
//! partial-order solver and the literal spec run on chain instances against this
//! module) and `partial_order::algorithm::test` (general posets).

use crate::structures::Direction;
use crate::total_order::stochastic_dominance::censored::structures::CensoredSdContext;

pub(crate) fn algorithm<D: Direction, X: crate::Float, Y: crate::Float>(
    context: &CensoredSdContext<X, Y>,
) -> Vec<f32> {
    (0..context.thresholds.len())
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
                            // Reference algorithm — accumulate in f64 for maximum precision;
                            // we downcast to f32 at the output to match the production cdfs.
                            let mut survival: f64 = 1.0;
                            for (i_rs, observation) in data.iter().enumerate() {
                                if threshold < observation.y {
                                    break;
                                }

                                if observation.observed {
                                    survival *= 1.0
                                        - observation.weight as f64
                                            / data[i_rs..]
                                                .iter()
                                                .map(|o| o.weight as f64)
                                                .sum::<f64>();
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
            // Clip in order of increasing span, so that every value used for clipping —
            // the same-row left part (r..=k) as well as the column right part (k+1..=s),
            // both spanning strictly fewer covariates — is itself already fully clipped.
            for span in 1..context.n_covariate() {
                for r in 0..context.n_covariate() - span {
                    let s = r + span;

                    // The clip intervals of the different splits k are always compatible:
                    // L = max_k min(left, right) never exceeds U = min_k max(left, right).
                    // Consequently, sequentially clipping into each split's interval (below)
                    // is equivalent to a single clamp into [L, U], independent of the k
                    // order — which is what the fast algorithm's bound propagation relies
                    // on.
                    let mut lower = f64::NEG_INFINITY;
                    let mut upper = f64::INFINITY;
                    for k in r..s {
                        let left = survivals[r][k - r];
                        let right = survivals[k + 1][s - k - 1];
                        lower = lower.max(left.min(right));
                        upper = upper.min(left.max(right));
                    }
                    assert!(
                        lower <= upper,
                        "conflicting clip bounds for interval ({r}, {s}): L={lower} > U={upper}",
                    );

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
        .map(|s| (1.0 - s) as f32)
        .collect()
}
