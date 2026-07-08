//! Large-instance tests for the censored stochastic-dominance fast path (`super::fast`) that drive
//! its SIMD kernels and acceleration machinery through the regimes small inputs never reach.
//!
//! The differential suites in [`super::test`] establish correctness, but only up to a few dozen
//! covariates. At that scale none of the *acceleration-relevant regimes* actually run: the wide
//! AVX-512 reduction loops (which only engage once a `propagate_bounds` row is >= 32 elements), the
//! 64-wide K-M `walk_scan` round, the cross-block pooling of hundreds of covariates, the
//! bound-propagation cache tiling (`POOL_TILE = 96`), and the pair-scheduled cell sweep. These
//! tests deliberately use **many unique covariates and very few thresholds**, so the estimator
//! triangle (`O(C^2)` cells, spans up to `C`) is large while the *answer* stays small and, crucially,
//! has a **hand-derived closed form** — the reference `super::definition` is `O(C^2 * n + C^3)` and
//! becomes the bottleneck well before `fast` does, so it cannot be the oracle at this scale.
//!
//! Two families, each isolating a different part of the machinery:
//!
//! * [`single_global_block`] — fully anti-monotone data collapses the entire grid into ONE pooled
//!   block whose value is the global Kaplan-Meier survival, so the CDF is constant across
//!   covariates with a closed form `(C+1)/(C+7)`, `(C+4)/(C+7)`. This drives the widest possible
//!   spans: the 32-wide AVX-512 `propagate_bounds` main loop, the 64-wide `walk_scan`, the
//!   pair-scheduled sweep, the `penultimate_start > split_x` skip, the uncensored→generalized bridge
//!   and the exact-0 completion pinning — all at hundreds of covariates.
//!
//! * [`wide_ultimate_merge`] — a narrow left region merges with a wide (> `POOL_TILE`) right block
//!   in a single pooling step, the one regime the anti-monotone collapse never reaches (its
//!   incremental merges always have a width-1 *ultimate*). This is what exercises the second tile of
//!   the `pool` cache-tiling loop. Verified against `super::definition` (still feasible at this
//!   width) rather than a closed form.

use crate::structures::Increasing;
use crate::total_order::preprocessing::preprocess_censored as preprocess;
use crate::total_order::stochastic_dominance::censored::{definition, fast};

/// Largest tolerated `|fast - reference|`. The fast path accumulates the Kaplan-Meier product and
/// the clip bounds in f32 with a different operation order than the reference, so agreement is only
/// expected to a few f32 ulps — the same bound the differential suite uses.
const TOL: f32 = 1e-5;

/// Build Construction A at `c` covariates (see [`single_global_block`] for the derivation).
///
/// Covariate `i` (its value is `i`, so the covariates are already sorted and unique) carries four
/// observations, all sharing the same tiny threshold grid `{y=1, y=2}` (the `y=3` survivors are
/// censored above the last event and create no threshold):
///   event    @ y=1  (w = i+1)   censored @ y=1  (w = 1)
///   event    @ y=2  (w = 1)     survivor @ y=3  (w = 1)
fn anti_monotone(c: usize) -> (Vec<f64>, Vec<f64>, Vec<bool>, Vec<f64>) {
    let mut x = Vec::with_capacity(4 * c);
    let mut y = Vec::with_capacity(4 * c);
    let mut observed = Vec::with_capacity(4 * c);
    let mut weight = Vec::with_capacity(4 * c);
    let mut push = |xi: usize, yi: f64, obs: bool, wi: f64| {
        x.push(xi as f64);
        y.push(yi);
        observed.push(obs);
        weight.push(wi);
    };
    for i in 0..c {
        push(i, 1.0, true, (i + 1) as f64); // event @ y=1, mass grows with i
        push(i, 1.0, false, 1.0); // censored @ y=1 (drives the generalized phase)
        push(i, 2.0, true, 1.0); // event @ y=2
        push(i, 3.0, false, 1.0); // survivor, censored above the last event
    }
    (x, y, observed, weight)
}

/// Construction A: strictly stochastically decreasing per-covariate sub-distributions force the
/// increasing S-IDR isotonic constraint to be violated between *every* pair of adjacent covariates,
/// so the pool-adjacent-violators solution collapses the whole grid into a **single block**. The
/// block value is the global Kaplan-Meier survival over all observations, so the CDF is constant
/// across covariates, with a closed form.
///
/// Per-covariate survival is `S_i(y1) = 3/(i+4)` and `S_i(y2) = 1.5/(i+4)`, both strictly
/// decreasing in `i` — hence the full pool. Because the events at a shared response are tied, their
/// Kaplan-Meier factors telescope (`prod (1 - w_k / R_k) = 1 - (sum w_k)/R`), so the *global*
/// survival is a two-factor product over the pooled risk sets:
///
/// ```text
///   R1 = sum_{i<C} (i+4) = C(C+7)/2     total weight = risk at y1
///   E1 = sum_{i<C} (i+1) = C(C+1)/2     event mass at y1        => S(y1) = 1 - E1/R1 = 6/(C+7)
///   R2 = R1 - E1 - C = 2C               risk at y2 (y1 events + y1 censoring have left)
///   E2 = C                              event mass at y2        => S(y2) = S(y1)(1 - E2/R2) = 3/(C+7)
/// ```
///
/// giving, at **every** covariate,
///
/// ```text
///   CDF(y1) = 1 - 6/(C+7) = (C+1)/(C+7)
///   CDF(y2) = 1 - 3/(C+7) = (C+4)/(C+7)
/// ```
///
/// Small sizes are cross-checked against the executable spec `definition`; the large size (well past
/// where the `O(C^2 n)` reference is comfortable, and deep into the multi-hundred-span SIMD regime)
/// is checked against the closed form alone. The same closed form has been verified against `fast`
/// up to `C = 1500` during development — the formula is exact at any size; only `fast`'s debug-build
/// runtime keeps the committed instance moderate.
#[test]
fn single_global_block() {
    // (C, cross-check against the O(C^3) reference?)  The reference is skipped at the large size.
    let cases = [(8usize, true), (24, true), (48, true), (128, false)];
    for (c, check_reference) in cases {
        let (x, y, observed, weight) = anti_monotone(c);
        let ctx = preprocess(&x, &y, &observed, &weight).unwrap();
        assert_eq!(ctx.n_covariate(), c);
        assert_eq!(ctx.thresholds, vec![1.0, 2.0], "C={c}");

        let fast_cdf = fast::algorithm::<Increasing, _, _>(&ctx, &crate::NoProgress);
        assert_eq!(fast_cdf.len(), 2 * c);

        // Closed form: constant across covariates, one value per threshold.
        let cdf_y1 = (c as f64 + 1.0) / (c as f64 + 7.0);
        let cdf_y2 = (c as f64 + 4.0) / (c as f64 + 7.0);
        for (t, expected) in [cdf_y1, cdf_y2].into_iter().enumerate() {
            let row = &fast_cdf[t * c..(t + 1) * c];
            for (i, &got) in row.iter().enumerate() {
                assert!(
                    (got as f64 - expected).abs() < TOL as f64,
                    "C={c} threshold={t} covariate={i}: got {got}, want {expected} (single-block CDF)",
                );
            }
        }

        // The whole grid pooled: rows must be exactly flat (single-block signature).
        for t in 0..2 {
            let row = &fast_cdf[t * c..(t + 1) * c];
            let (lo, hi) = row.iter().fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
            assert!(hi - lo < TOL, "C={c} threshold={t}: row not constant, spread {}", hi - lo);
        }
        // Proper CDF: nondecreasing along thresholds at every covariate.
        assert!(cdf_y1 <= cdf_y2);

        if check_reference {
            let def_cdf = definition::algorithm::<Increasing, _, _>(&ctx);
            let worst = fast_cdf
                .iter()
                .zip(&def_cdf)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f32, f32::max);
            assert!(worst < TOL, "C={c}: fast vs definition worst {worst:e}");
        }
    }
}

/// Build Construction B at right-block width `w` (see [`wide_ultimate_merge`]).
///
/// Two narrow left singletons `{0, 1}` with moderate `y=1` event mass, then a wide right plateau
/// `{2 ..= w+1}` whose covariates have tiny, strictly increasing `y=1` event mass (so the plateau
/// pools into ONE block at the uncensored `y=1` threshold) plus heavy `y=1` censoring. The
/// right-most plateau covariate additionally carries a huge `y=2` event whose weight scales with
/// `w`, so once it is applied the plateau's survival crashes below the left singletons.
fn wide_plateau(w: usize) -> (Vec<f64>, Vec<f64>, Vec<bool>, Vec<f64>) {
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut observed = Vec::new();
    let mut weight = Vec::new();
    let mut push = |xi: usize, yi: f64, obs: bool, wi: f64| {
        x.push(xi as f64);
        y.push(yi);
        observed.push(obs);
        weight.push(wi);
    };

    // Left narrow region: two singletons with moderate, decreasing-in-value event mass. Their
    // survival stays below the plateau's at y=1 (no early violation) but above it after the crash.
    for (k, ev) in [(0usize, 5.0), (1, 4.0)] {
        push(k, 1.0, true, ev); // event @ y=1
        push(k, 2.0, true, 1.0); // event @ y=2 (small)
        push(k, 3.0, false, 1.0); // survivor
    }

    // Wide right plateau.
    let plateau_end = 2 + w; // exclusive covariate index
    let crash_covariate = plateau_end - 1; // right-most plateau covariate
    for cov in 2..plateau_end {
        let j = cov - 2;
        // Tiny, strictly increasing event mass -> survival slightly decreasing -> the plateau pools.
        push(cov, 1.0, true, 1.0 + j as f64 * 0.02);
        push(cov, 1.0, false, 5.0); // heavy censoring at y=1
        push(cov, 3.0, false, 1.0); // survivor
        if cov == crash_covariate {
            // Huge y=2 event: crashes the pooled plateau survival below the left singletons.
            push(cov, 2.0, true, 10.0 * w as f64);
        }
    }
    (x, y, observed, weight)
}

/// Construction B: forces a single `pool` merge whose *ultimate* (right / column) block spans more
/// than `POOL_TILE` (= 96) covariates, driving the bound-propagation cache-tiling loop through more
/// than one tile — the regime [`single_global_block`] never reaches (its incremental collapse only
/// ever merges width-1 ultimates).
///
/// The wide right plateau pools into one block at the uncensored `y=1` threshold; entering the
/// generalized phase it is a single wide partition. The right-most plateau covariate's huge `y=2`
/// event then drops the plateau's survival below the left singletons, so the whole plateau — as the
/// wide *ultimate* of the merge — is pooled left in one step over `w + 1` columns (> 96), spanning
/// several `POOL_TILE`-sized tiles.
///
/// No closed form is used: the answer is checked against the executable specification
/// `super::definition`, exactly like the differential suites, but on an instance hand-built to reach
/// the multi-tile path. `w = 140` keeps the wide merge (141 columns) comfortably above `POOL_TILE`
/// even if that constant is retuned within its documented 64–128 range.
#[test]
fn wide_ultimate_merge() {
    let w = 140usize;
    let (x, y, observed, weight) = wide_plateau(w);
    let ctx = preprocess(&x, &y, &observed, &weight).unwrap();
    assert_eq!(ctx.n_covariate(), w + 2);
    assert_eq!(ctx.thresholds, vec![1.0, 2.0]);

    let fast_cdf = fast::algorithm::<Increasing, _, _>(&ctx, &crate::NoProgress);
    let def_cdf = definition::algorithm::<Increasing, _, _>(&ctx);
    assert_eq!(fast_cdf.len(), def_cdf.len());

    let worst = fast_cdf
        .iter()
        .zip(&def_cdf)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f32, f32::max);
    assert!(worst < TOL, "fast vs definition worst {worst:e}");

    let nc = ctx.n_covariate();
    // Structural signature we engineered: at y=1 the plateau is one flat block; at y=2 the crash has
    // pooled the entire grid into one flat block.
    let y1 = &fast_cdf[..nc];
    let plateau = &y1[2..nc - 1]; // covariates 2 ..= w  (the pooled plateau, excl. the crash cov)
    let (plo, phi) = plateau
        .iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &v| (lo.min(v), hi.max(v)));
    assert!(phi - plo < TOL, "y=1 plateau not flat: spread {}", phi - plo);

    let y2 = &fast_cdf[nc..2 * nc];
    let (lo, hi) = y2
        .iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &v| (lo.min(v), hi.max(v)));
    assert!(hi - lo < TOL, "y=2 row not fully pooled: spread {}", hi - lo);

    // Every threshold row of an increasing fit is nonincreasing along the sorted covariates.
    for (t, row) in fast_cdf.chunks_exact(nc).enumerate() {
        assert!(
            row.iter().rev().is_sorted(),
            "threshold {t}: row must be nonincreasing along covariates",
        );
    }
}
