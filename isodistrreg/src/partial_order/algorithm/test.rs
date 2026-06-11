//! Anchor and differential tests for the censored partial-order solver against the
//! executable specification in [`super::definition`].
//!
//! The production solver reduces the per-threshold problem to greedy extreme-lower-set
//! peeling (`partial_order::functionals::algorithm_pre_sorted_inner`); the spec computes
//! the defining min-max formula literally and asserts the min-max/max-min saddle
//! equality internally. Agreement here is the evidence that the greedy construction
//! solves the formula.

use crate::NoProgress;
use crate::partial_order::{self, CensoredContext, CovariateGroups};
use crate::structures::{Decreasing, Increasing};

/// Same comparison tolerance as the total-order censored differential suite: both
/// pipelines narrow to f32 (the solver computes its clipping triangle in f32, the spec
/// narrows once at output), so disagreements are bounded by f32 accumulation noise.
const TOLERANCE: f32 = 1e-5;

type Instance = (usize, Vec<f64>, Vec<f64>, Vec<bool>, Vec<f64>);

fn preprocess(instance: &Instance) -> CensoredContext<f64, f64> {
    let (dimension, x, y, observed, weights) = instance;
    partial_order::preprocess_censored(x, y, observed, weights, &CovariateGroups::empty(*dimension))
}

fn compare(label: &str, actual: &[f32], expected: &[f32]) -> Result<(), String> {
    if actual.len() != expected.len() {
        return Err(format!(
            "{label}: length mismatch — solver {} vs spec {}",
            actual.len(),
            expected.len()
        ));
    }
    let worst = actual
        .iter()
        .zip(expected)
        .map(|(a, e)| (a - e).abs())
        .fold(0.0f32, f32::max);
    if worst <= TOLERANCE {
        Ok(())
    } else {
        Err(format!(
            "{label}: max diff {worst:e}\n  solver: {actual:?}\n  spec:   {expected:?}"
        ))
    }
}

/// Structural sanity on any output: every threshold row respects the covariate order
/// (Increasing: comparable i ⪯ j implies F_i >= F_j) and every covariate column is
/// nondecreasing across thresholds. A small epsilon absorbs f32 noise; genuine
/// violations were audit findings, not noise.
fn assert_structure(
    cdfs: &[f32],
    context: &CensoredContext<f64, f64>,
    increasing: bool,
    label: &str,
) {
    let m = context.n_covariate();
    if m == 0 {
        return;
    }
    let dimension = context.covariate_dimension();
    let leq = |i: usize, j: usize| {
        (0..dimension)
            .all(|k| context.x_unique[i * dimension + k] <= context.x_unique[j * dimension + k])
    };
    for (t, row) in cdfs.chunks_exact(m).enumerate() {
        for i in 0..m {
            for j in 0..m {
                if i != j && leq(i, j) {
                    let (hi, lo) = if increasing { (i, j) } else { (j, i) };
                    assert!(
                        row[hi] >= row[lo] - 1e-6,
                        "{label}: covariate order violated at threshold {t}: \
                         F[{hi}] = {} < F[{lo}] = {}",
                        row[hi],
                        row[lo],
                    );
                }
            }
        }
    }
    for i in 0..m {
        for t in 1..cdfs.len() / m {
            assert!(
                cdfs[t * m + i] >= cdfs[(t - 1) * m + i] - 1e-6,
                "{label}: column {i} decreases between thresholds {} and {t}",
                t - 1,
            );
        }
    }
}

fn check_both_directions(instance: &Instance) -> Result<(), String> {
    let context = preprocess(instance);

    let solver_inc = partial_order::censored::<Increasing, _, _>(&context, &NoProgress);
    let spec_inc = partial_order::censored_definition::<Increasing, _, _>(&context);
    compare("Increasing", &solver_inc.cdfs, &spec_inc)?;
    assert_structure(&solver_inc.cdfs, &context, true, "solver Increasing");
    assert_structure(&spec_inc, &context, true, "spec Increasing");

    let solver_dec = partial_order::censored::<Decreasing, _, _>(&context, &NoProgress);
    let spec_dec = partial_order::censored_definition::<Decreasing, _, _>(&context);
    compare("Decreasing", &solver_dec.cdfs, &spec_dec)?;
    assert_structure(&solver_dec.cdfs, &context, false, "solver Decreasing");
    assert_structure(&spec_dec, &context, false, "spec Decreasing");

    Ok(())
}

/// Incomparable covariates impose no constraints: each keeps its own Kaplan-Meier
/// curve. Covariate A = (1, 2): event at y = 1, censored at y = 3 (so its curve stays
/// at 1/2 without completing); B = (2, 1): event at y = 2.
#[test]
fn antichain_keeps_own_kaplan_meier() {
    let instance: Instance = (
        2,
        vec![1.0, 2.0, 1.0, 2.0, 2.0, 1.0],
        vec![1.0, 3.0, 2.0],
        vec![true, false, true],
        vec![1.0, 1.0, 1.0],
    );
    let context = preprocess(&instance);
    let expected = [0.5, 0.0, 0.5, 1.0];

    let solver = partial_order::censored::<Increasing, _, _>(&context, &NoProgress);
    assert_eq!(solver.cdfs, expected);
    let spec = partial_order::censored_definition::<Increasing, _, _>(&context);
    assert_eq!(spec, expected);
}

/// A comparable pair with antitone data pools (Increasing fit): A = (1,1) ⪯ B = (2,2),
/// A's event at y = 2, B's at y = 1, so at the first threshold the empirical values
/// (0, 1) violate F_A >= F_B and both take RKM_{A,B} = 1/2 = the pooled Kaplan-Meier;
/// at the second threshold both curves complete at exactly 1.
#[test]
fn comparable_pair_pools_under_increasing() {
    let instance: Instance = (
        2,
        vec![1.0, 1.0, 2.0, 2.0],
        vec![2.0, 1.0],
        vec![true, true],
        vec![1.0, 1.0],
    );
    let context = preprocess(&instance);

    let solver = partial_order::censored::<Increasing, _, _>(&context, &NoProgress);
    assert_eq!(solver.cdfs, [0.5, 0.5, 1.0, 1.0]);
    let spec = partial_order::censored_definition::<Increasing, _, _>(&context);
    assert_eq!(spec, [0.5, 0.5, 1.0, 1.0]);
}

/// The same data is feasible for a Decreasing fit (F isotone along the order), so
/// both covariates keep their own Kaplan-Meier values.
#[test]
fn comparable_pair_is_feasible_under_decreasing() {
    let instance: Instance = (
        2,
        vec![1.0, 1.0, 2.0, 2.0],
        vec![2.0, 1.0],
        vec![true, true],
        vec![1.0, 1.0],
    );
    let context = preprocess(&instance);

    let solver = partial_order::censored::<Decreasing, _, _>(&context, &NoProgress);
    assert_eq!(solver.cdfs, [0.0, 1.0, 1.0, 1.0]);
    let spec = partial_order::censored_definition::<Decreasing, _, _>(&context);
    assert_eq!(spec, [0.0, 1.0, 1.0, 1.0]);
}

mod differential {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    /// dim-D integer-grid covariates under the componentwise order: small grids
    /// maximize duplicate rows, incomparabilities, response ties, and censoring
    /// interactions.
    fn integer_grid(
        rng: &mut StdRng,
        n: usize,
        dimension: usize,
        covariate_levels: u32,
        response_levels: u32,
        censored_percent: u32,
        random_weights: bool,
    ) -> Instance {
        let x = (0..n * dimension)
            .map(|_| rng.random_range(1..=covariate_levels) as f64)
            .collect();
        let y = (0..n)
            .map(|_| rng.random_range(1..=response_levels) as f64)
            .collect();
        let mut observed: Vec<bool> = (0..n)
            .map(|_| rng.random_range(0..100) >= censored_percent)
            .collect();
        // The uncensored case is covered by the total-order suites; force at least
        // one censored observation so every instance exercises the censored paths.
        if observed.iter().all(|&o| o) {
            let flip = rng.random_range(0..n);
            observed[flip] = false;
        }
        let weights = (0..n)
            .map(|_| {
                if random_weights {
                    rng.random_range(0.5..2.0)
                } else {
                    1.0
                }
            })
            .collect();
        (dimension, x, y, observed, weights)
    }

    /// The spec's outer reduction enumerates all (lower set, upper set) pairs; cap
    /// that product so adversarially antichain-like draws cannot blow up the debug
    /// runtime. The guard mirrors the spec's own componentwise order construction.
    fn within_budget(context: &CensoredContext<f64, f64>) -> bool {
        let m = context.n_covariate();
        if m == 0 {
            // All-censored draw: the empty context is trivially comparable.
            return true;
        }
        if m > 12 {
            return false;
        }
        let dimension = context.covariate_dimension();
        let leq = |i: usize, j: usize| {
            (0..dimension)
                .all(|k| context.x_unique[i * dimension + k] <= context.x_unique[j * dimension + k])
        };
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
        let full: u64 = (1 << m) - 1;
        let mut lowers = 0usize;
        let mut uppers = 0usize;
        for set in 1..=full {
            if (0..m).all(|j| set & (1 << j) == 0 || below[j] & !set == 0) {
                lowers += 1;
            }
            if (0..m).all(|j| set & (1 << j) == 0 || above[j] & !set == 0) {
                uppers += 1;
            }
        }
        lowers * uppers <= 30_000
    }

    /// Sweep small random posets through both the production solver and the literal
    /// specification, in both directions. The saddle equality is asserted inside the
    /// spec on every evaluation; the seeded generator makes any panic reproducible.
    #[test]
    fn poset_sweep() {
        let mut failures = Vec::new();
        let mut ran = 0usize;
        let mut skipped = 0usize;

        let mut seed = 0u64;
        for dimension in [2usize, 3] {
            for covariate_levels in [2u32, 3, 4] {
                for response_levels in [2u32, 3, 5] {
                    for censored_percent in [20u32, 40, 70] {
                        for n in [4usize, 6, 8, 10, 12, 14] {
                            for _ in 0..10 {
                                seed += 1;
                                let mut rng = StdRng::seed_from_u64(seed);
                                let instance = integer_grid(
                                    &mut rng,
                                    n,
                                    dimension,
                                    covariate_levels,
                                    response_levels,
                                    censored_percent,
                                    seed.is_multiple_of(2),
                                );
                                let context = preprocess(&instance);
                                if !within_budget(&context) {
                                    skipped += 1;
                                    continue;
                                }
                                ran += 1;
                                if let Err(message) = check_both_directions(&instance) {
                                    failures.push(format!(
                                        "seed {seed} (dim {dimension}, n {n}, levels \
                                         {covariate_levels}/{response_levels}, censoring \
                                         {censored_percent}%): {message}\n  instance: {instance:?}"
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }

        assert!(
            failures.is_empty(),
            "{} of {ran} instances diverged:\n{}",
            failures.len(),
            failures.join("\n"),
        );
        // The budget guard must stay an edge-case escape hatch, not silently hollow
        // out the suite.
        assert!(
            ran * 10 >= (ran + skipped) * 9,
            "budget guard skipped too many instances: {skipped} of {}",
            ran + skipped,
        );
    }
}
