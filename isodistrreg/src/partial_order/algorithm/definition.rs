//! Executable specification of the censored partial-order IDR estimator.
//!
//! Implements the defining formula literally, with no algorithmic shortcuts: for
//! unique covariates ξ_1, .., ξ_m under a partial order with lower sets 𝓛 and upper
//! sets 𝓤,
//!
//! ```text
//! F̂_{ξ_i}(y) = min_{U ∈ 𝓤: i ∈ U}  max_{L ∈ 𝓛: i ∈ L}  RKM_{L ∩ U}(y)
//! RKM_S(y)   = clamp( KM_S(y),
//!                     max_{(L,U) ∈ P_S} min(RKM_{S∩L}, RKM_{S∩U}),
//!                     min_{(L,U) ∈ P_S} max(RKM_{S∩L}, RKM_{S∩U}) )
//! P_S        = { (S∩L, S∩U) : L ∈ 𝓛, U ∈ 𝓤, L ∪ U ⊇ S, L ∩ U = ∅ }
//! ```
//!
//! where KM_S is the weighted Kaplan-Meier estimator over the observations of the
//! covariates in S (events before censorings at ties), and RKM is the recursively
//! clipped self-consistent Kaplan-Meier estimator. Degenerate pairs (S, ∅)/(∅, S)
//! are excluded from P_S — RKM_∅ is undefined — matching the production convention
//! (`functionals::recursive_clipping` skips the empty/full split, and the total-order
//! definition only considers interior split points).
//!
//! The implementation enumerates subsets as `u64` masks (so `m <= 64`; the suites cap
//! `m <= 12`) and derives the order relation directly from `x_unique` componentwise —
//! deliberately independent of `compute_transitive_closure`/`derive_transitive_reduction`,
//! so the differential suites cover those too. All arithmetic is f64, narrowed to f32
//! once per output value like the total-order definition.
//!
//! The outer reduction computes BOTH `min_U max_L` (the specification) and its dual
//! `max_L min_U` and asserts they are exactly equal: weak duality guarantees
//! `max-min <= min-max`, both sides select from the same finite table of f64 RKM
//! values with no further arithmetic, and the production greedy peeling
//! (`partial_order::functionals::algorithm_pre_sorted_inner`) is only equivalent to
//! the formula when this saddle equality holds. A gap is a mathematical finding and
//! must surface loudly, never be resolved silently.

use crate::partial_order::CensoredContext;
use crate::structures::Direction;
use std::collections::HashMap;

/// Brute-force reference: one CDF value per (threshold, unique covariate), threshold-major,
/// covariates in `x_unique` order — directly comparable to
/// [`partial_order::censored`](crate::partial_order::censored)'s `AlgorithmOutput::cdfs`.
pub(crate) fn algorithm<D: Direction, X: PartialOrd, Y>(
    context: &CensoredContext<X, Y>,
) -> Vec<f32> {
    if context.n() == 0 {
        return Vec::with_capacity(0);
    }

    let m = context.n_covariate();
    assert!(m <= 64, "the executable spec uses u64 subset masks");
    let dimension = context.covariate_dimension();

    // Componentwise order on the unique covariate rows. `x_unique` has already been
    // transformed by `unify_group_orders`, so the componentwise comparison is the
    // partial order the production solver uses.
    let leq = |i: usize, j: usize| {
        (0..dimension).all(|k| {
            let l = &context.x_unique[i * dimension + k];
            let r = &context.x_unique[j * dimension + k];
            l <= r
        })
    };

    // Inclusive closure masks: below[j] = { i : ξ_i ⪯ ξ_j }, above[j] = { i : ξ_j ⪯ ξ_i }.
    let mut below = vec![0u64; m];
    let mut above = vec![0u64; m];
    for j in 0..m {
        for i in 0..m {
            if leq(i, j) {
                below[j] |= 1 << i;
            }
            if leq(j, i) {
                above[j] |= 1 << i;
            }
        }
    }

    // Group observations by unique covariate: (threshold index, observed, weight).
    let mut observations = vec![Vec::new(); m];
    for j in 0..context.n() {
        observations[context.x[j]].push((context.y[j], context.y_observed[j], context.weight[j]));
    }

    // All global lower and upper sets, and per covariate the ones containing it.
    let full: u64 = if m == 64 { !0 } else { (1 << m) - 1 };
    let mut lowers = Vec::new();
    let mut uppers = Vec::new();
    for set in 1..=full {
        if (0..m).all(|j| set & (1 << j) == 0 || below[j] & !set == 0) {
            lowers.push(set);
        }
        if (0..m).all(|j| set & (1 << j) == 0 || above[j] & !set == 0) {
            uppers.push(set);
        }
    }

    let spec = Spec {
        below,
        above,
        observations,
    };

    // For a decreasing fit the covariate order is reversed, which swaps the lower-set
    // and upper-set families; RKM and P_S are direction-free (S ∩ L and S ∩ U swap
    // roles, but the clamp bounds are symmetric in the pair).
    let (min_family, max_family) = if D::IS_INCREASING {
        (&uppers, &lowers)
    } else {
        (&lowers, &uppers)
    };

    let mut cdfs = Vec::with_capacity(context.n_threshold() * m);
    for t in 0..context.n_threshold() {
        let mut memo = HashMap::new();
        for i in 0..m {
            let bit = 1u64 << i;

            // Primary orientation: the specification's min over the outer family of
            // the max over the inner family. `i ∈ L ∩ U` keeps every intersection
            // nonempty.
            let primary = min_family
                .iter()
                .filter(|&&outer| outer & bit != 0)
                .map(|&outer| {
                    max_family
                        .iter()
                        .filter(|&&inner| inner & bit != 0)
                        .map(|&inner| spec.rkm(outer & inner, t, &mut memo))
                        .fold(f64::NEG_INFINITY, f64::max)
                })
                .fold(f64::INFINITY, f64::min);

            // Dual orientation (max over inner family of min over outer family). Weak
            // duality gives dual <= primary; both select from the same finite RKM
            // table, so exact equality is the saddle condition the greedy production
            // solver relies on.
            let dual = max_family
                .iter()
                .filter(|&&inner| inner & bit != 0)
                .map(|&inner| {
                    min_family
                        .iter()
                        .filter(|&&outer| outer & bit != 0)
                        .map(|&outer| spec.rkm(outer & inner, t, &mut memo))
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(f64::NEG_INFINITY, f64::max);

            assert!(
                primary == dual,
                "saddle gap in the min-max formula at covariate {i}, threshold {t}: \
                 min-max = {primary}, max-min = {dual} — the greedy poset solver is \
                 not equivalent to the specification on this instance",
            );

            cdfs.push(primary as f32);
        }
    }
    cdfs
}

struct Spec {
    below: Vec<u64>,
    above: Vec<u64>,
    /// Per unique covariate: (threshold index, observed, weight).
    observations: Vec<Vec<(usize, bool, f64)>>,
}

impl Spec {
    /// Inclusive downward closure of a subset.
    fn down(&self, set: u64) -> u64 {
        let mut result = 0;
        let mut rest = set;
        while rest != 0 {
            let j = rest.trailing_zeros() as usize;
            result |= self.below[j];
            rest &= rest - 1;
        }
        result
    }

    /// Inclusive upward closure of a subset.
    fn up(&self, set: u64) -> u64 {
        let mut result = 0;
        let mut rest = set;
        while rest != 0 {
            let j = rest.trailing_zeros() as usize;
            result |= self.above[j];
            rest &= rest - 1;
        }
        result
    }

    /// Weighted Kaplan-Meier CDF of the covariates in `set` at threshold index `t`,
    /// mirroring `functionals::KaplanMeier::evaluate` in f64: every observation of the
    /// set is at risk, observations with response at or below the threshold leave the
    /// risk set in (response, events-before-censorings) order, events multiply the
    /// survival, and a survival that is exactly 0 in exact arithmetic (all mass at or
    /// below the threshold, last observation an event) is pinned.
    fn km(&self, set: u64, t: usize) -> f64 {
        let mut at_risk = 0.0;
        let mut n_total = 0usize;
        let mut subset = Vec::new();
        let mut rest = set;
        while rest != 0 {
            let j = rest.trailing_zeros() as usize;
            for &(y, observed, weight) in &self.observations[j] {
                at_risk += weight;
                n_total += 1;
                if y <= t {
                    subset.push((y, observed, weight));
                }
            }
            rest &= rest - 1;
        }

        subset.sort_unstable_by(|l, r| l.0.cmp(&r.0).then(l.1.cmp(&r.1).reverse()));

        let completes = subset.len() == n_total && subset.last().is_some_and(|o| o.1);

        let mut survival = 1.0;
        for (_, observed, weight) in subset {
            if observed {
                survival *= 1.0 - weight / at_risk;
            }
            at_risk -= weight;
        }
        if completes {
            survival = 0.0;
        }
        1.0 - survival
    }

    /// The recursively clipped Kaplan-Meier value of `set` at threshold index `t`.
    ///
    /// Enumerates every unordered two-block partition (a, b) of `set` and admits it as
    /// an element of P_S via the canonical-witness lemma: a global pair (L, U) with
    /// S ∩ L = a, S ∩ U = b, L ∩ U = ∅, L ∪ U ⊇ S exists if and only if
    /// down(a) ∩ up(b) = ∅. (⇐: L₀ = down(a) and U₀ = up(b) are a lower and an upper
    /// set with a ⊆ L₀, b ⊆ U₀, S ⊆ a ∪ b ⊆ L₀ ∪ U₀, and S ∩ L₀ = a / S ∩ U₀ = b
    /// because down(a) ∩ b ⊆ down(a) ∩ up(b) = ∅. ⇒: any valid L contains down(a) and
    /// any valid U contains up(b), so L ∩ U = ∅ forces down(a) ∩ up(b) = ∅.)
    ///
    /// The clamp bounds are symmetric in the pair's two values, so each unordered pair
    /// contributes once if either orientation admits a witness.
    fn rkm(&self, set: u64, t: usize, memo: &mut HashMap<u64, f64>) -> f64 {
        debug_assert_ne!(set, 0, "RKM of the empty set is undefined");
        if let Some(&value) = memo.get(&set) {
            return value;
        }

        let raw = self.km(set, t);

        let mut lower = f64::NEG_INFINITY;
        let mut upper = f64::INFINITY;
        // Standard submask enumeration; visit each unordered proper pair once.
        let mut a = (set - 1) & set;
        while a != 0 {
            let b = set & !a;
            if a < b && (self.down(a) & self.up(b) == 0 || self.down(b) & self.up(a) == 0) {
                let value_a = self.rkm(a, t, memo);
                let value_b = self.rkm(b, t, memo);
                lower = lower.max(value_a.min(value_b));
                upper = upper.min(value_a.max(value_b));
            }
            a = (a - 1) & set;
        }

        let value = raw.max(lower).min(upper);
        memo.insert(set, value);
        value
    }
}
