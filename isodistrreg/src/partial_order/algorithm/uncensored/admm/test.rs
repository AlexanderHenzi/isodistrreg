//! Correctness tests for the in-tree solver.
//!
//! Every expected value here is either derivable by hand or produced by an implementation
//! that shares no code with the solver. The differential oracle is
//! [`crate::partial_order::tonic_regression_pre_sorted`], which computes the same
//! projection by greedy extreme-lower-set peeling -- a combinatorial algorithm with no
//! linear algebra in it at all, so agreement is evidence rather than a tautology.
//!
//! The bridge between the two is the objective identity
//!
//! ```text
//!     0.5 ||u||^2 + q^T u  =  0.5 ||u - (-q)||^2 - 0.5 ||q||^2
//! ```
//!
//! so the QP is exactly the Euclidean projection of `-q` onto the cone `{A u >= 0}`. When
//! the coefficients are `+-1` that cone is the isotone cone of the poset, and the answer is
//! the unweighted isotonic regression of `-q`.
use super::constraints::Constraints;
use super::hildreth::Hildreth;
use super::kkt::KktFactor;
use super::polish::{self, Polisher};
use super::solver::{self, Settings, Status};
use crate::functionals::Average;
use crate::partial_order::routines::derive_transitive_reduction;
use crate::partial_order::tonic_regression_pre_sorted;
use crate::structures::{Increasing, Observation};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

fn settings() -> Settings {
    Settings {
        verbose: false,
        eps_abs: 1e-10,
        eps_rel: 0.0,
        max_iter: 20_000,
    }
}

/// Run the solver plus polish on one instance, returning the primal solution.
fn solve(edges: &[(usize, usize)], coef: &[(f64, f64)], q: &[f64]) -> Vec<f64> {
    let n = q.len();
    let m = edges.len();
    let constraints = Constraints { edges, coef, n };
    let mut factor = KktFactor::new(n, edges);
    let mut workspace = solver::Workspace::new(n, m);
    let mut polisher = Polisher::new(n, m);

    let outcome = solver::solve(
        &constraints,
        q,
        &settings(),
        &mut workspace,
        &mut factor,
        &mut polisher,
    );
    assert_eq!(outcome.status, Status::Solved, "solver did not converge");

    let solution = if outcome.finished_exact {
        workspace.x.clone()
    } else {
        let report = polish::polish(
            &constraints,
            q,
            |r| workspace.clipped[r],
            &workspace.x,
            &mut polisher,
        );
        if report.accepted {
            polisher.u.clone()
        } else {
            workspace.x.clone()
        }
    };

    assert_certified_optimal(&constraints, q, &solution);
    solution
}

fn objective(u: &[f64], q: &[f64]) -> f64 {
    0.5 * u.iter().map(|v| v * v).sum::<f64>() + u.iter().zip(q).map(|(a, b)| a * b).sum::<f64>()
}

/// Assert primal feasibility and a certified duality gap.
///
/// The dual `max_{lambda >= 0} -0.5 ||A^T lambda - q||^2` lower-bounds the optimum at any
/// dual-feasible point, with no convergence assumption, so `f(u) - g(lambda)` is a genuine
/// bound on how suboptimal `u` is -- not a comparison against another approximation.
fn assert_certified_optimal(a: &Constraints<'_>, q: &[f64], u: &[f64]) {
    let mut row_values = vec![0.0; a.m()];
    a.mul(u, &mut row_values);
    let worst = row_values.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(worst >= -1e-9, "infeasible: min (A u) = {worst:e}");

    let mut certifier = Hildreth::new(a.n, a.m());
    certifier.seed_zero(a, q);
    certifier.run(a, 5_000);
    let gap = objective(u, q) - certifier.dual_bound();
    assert!(gap <= 1e-7, "certified gap {gap:e} too large");
}

/// With no constraints the minimizer of `0.5 u^T u + q^T u` is `u = -q`.
#[test]
fn antichain_is_the_unconstrained_minimizer() {
    let q = [0.75, -2.0, 0.5, 3.25];
    let solution = solve(&[], &[], &q);
    for (&got, &q_c) in solution.iter().zip(&q) {
        assert!((got + q_c).abs() < 1e-9, "{solution:?}");
    }
}

/// Two nodes, `u_1 >= u_0`, with `q = (-1, 1)` so the unconstrained answer `(1, -1)`
/// violates it. The constraint binds, both pool to `t` minimizing `t^2 + (q_0 + q_1) t`,
/// giving `t = -(q_0 + q_1)/2 = 0`.
#[test]
fn single_binding_edge_pools_to_the_midpoint() {
    let solution = solve(&[(0, 1)], &[(-1.0, 1.0)], &[-1.0, 1.0]);
    assert!(solution[0].abs() < 1e-9, "{solution:?}");
    assert!(solution[1].abs() < 1e-9, "{solution:?}");
}

/// Data already satisfying the order is its own projection, exactly.
#[test]
fn feasible_data_is_its_own_projection() {
    // -q = (0, 1, 2, 3) is isotone along the chain, so u = -q.
    let edges: Vec<(usize, usize)> = (0..3).map(|i| (i, i + 1)).collect();
    let coef = vec![(-1.0, 1.0); 3];
    let q = [0.0, -1.0, -2.0, -3.0];
    let solution = solve(&edges, &coef, &q);
    for (&got, &q_c) in solution.iter().zip(&q) {
        assert!((got + q_c).abs() < 1e-9, "{solution:?}");
    }
}

/// Strictly anti-monotone data on a chain pools everything to the global mean.
#[test]
fn anti_monotone_chain_pools_to_the_global_mean() {
    let n = 5;
    let edges: Vec<(usize, usize)> = (0..n - 1).map(|i| (i, i + 1)).collect();
    let coef = vec![(-1.0, 1.0); n - 1];
    // -q = (4, 3, 2, 1, 0): strictly decreasing along the chain.
    let q: Vec<f64> = (0..n).map(|i| i as f64 - 4.0).collect();
    let expected = -q.iter().sum::<f64>() / n as f64;
    let solution = solve(&edges, &coef, &q);
    for &got in &solution {
        assert!((got - expected).abs() < 1e-9, "{solution:?} vs {expected}");
    }
}

/// Random 2-D posets against the combinatorial oracle.
///
/// `n <= 8` keeps the oracle's ideal enumeration bounded by `2^8`, so the sweep stays
/// cheap enough for CI while still covering antichains, chains and everything between.
#[test]
fn matches_the_exact_poset_oracle() {
    let mut worst = 0.0f64;
    let mut instances = 0usize;
    for seed in 0..200u64 {
        let mut rng = StdRng::seed_from_u64(seed);
        let n_request = rng.random_range(2..=8usize);
        // Random 2-D covariates, deduplicated and lexicographically sorted, exactly the
        // shape `derive_transitive_reduction` is fed in production.
        let levels = rng.random_range(2..=4i64);
        let mut points: Vec<(i64, i64)> = (0..n_request)
            .map(|_| (rng.random_range(0..levels), rng.random_range(0..levels)))
            .collect();
        points.sort_unstable();
        points.dedup();
        let n = points.len();
        if n < 2 {
            continue;
        }
        let flat: Vec<f64> = points
            .iter()
            .flat_map(|&(a, b)| [a as f64, b as f64])
            .collect();
        let edges = derive_transitive_reduction(&flat, n, 2);
        let coef = vec![(-1.0, 1.0); edges.len()];
        let q: Vec<f64> = (0..n).map(|_| rng.random_range(-3.0..3.0f64)).collect();

        let solution = solve(&edges, &coef, &q);

        // The oracle projects `-q` onto the isotone cone of the same edge set.
        let targets: Vec<f64> = q.iter().map(|value| -value).collect();
        let expected = tonic_regression_pre_sorted::<Increasing, _, _>(
            targets.iter().copied(),
            &edges,
            &Average,
        );
        for (&got, &want) in solution.iter().zip(&expected) {
            worst = worst.max((got - want).abs());
        }
        assert!(
            worst < 1e-6,
            "seed {seed}: solver {solution:?} vs oracle {expected:?} (points {points:?})"
        );
        instances += 1;
    }
    assert!(instances >= 150, "only {instances} usable instances");
    println!("exact-oracle sweep: {instances} instances, worst deviation {worst:e}");
}

/// The same differential test with unequal weights, which is what makes the coefficients
/// differ from `+-1` -- the case the polish's ratio propagation has to get right and the
/// one real dumped instances never exercise.
#[test]
fn weighted_instances_match_the_weighted_oracle() {
    let mut worst = 0.0f64;
    for seed in 1000..1120u64 {
        let mut rng = StdRng::seed_from_u64(seed);
        let levels = rng.random_range(2..=4i64);
        let mut points: Vec<(i64, i64)> = (0..rng.random_range(2..=8usize))
            .map(|_| (rng.random_range(0..levels), rng.random_range(0..levels)))
            .collect();
        points.sort_unstable();
        points.dedup();
        let n = points.len();
        if n < 2 {
            continue;
        }
        let flat: Vec<f64> = points
            .iter()
            .flat_map(|&(a, b)| [a as f64, b as f64])
            .collect();
        let edges = derive_transitive_reduction(&flat, n, 2);

        // Weights spanning several orders of magnitude, as powers of two so the
        // energy-coordinate rescaling is exact.
        let weights: Vec<f64> = (0..n)
            .map(|_| 2f64.powi(rng.random_range(-4..4i32)))
            .collect();
        let targets: Vec<f64> = (0..n).map(|_| rng.random_range(0.0..1.0f64)).collect();

        // Energy coordinates, mirroring `algorithm` exactly -- including the power-of-two
        // `weight_scale` that places the largest weight in (1/2, 1], which is what keeps
        // every sqrt_weight at or below 1. Omitting it produces constraint rows the
        // production path never generates.
        let max_weight = weights.iter().copied().fold(0.0f64, f64::max);
        let mut weight_scale = 2f64.powi(-(max_weight.log2().ceil() as i32).clamp(-1022, 1022));
        if max_weight * weight_scale > 1.0 {
            weight_scale *= 0.5;
        }
        let sqrt_weight: Vec<f64> = weights.iter().map(|w| (w * weight_scale).sqrt()).collect();
        let q_energy: Vec<f64> = weights
            .iter()
            .zip(&targets)
            .zip(&sqrt_weight)
            .map(|((w, p), sw)| -weight_scale * w * p / sw)
            .collect();
        let coef: Vec<(f64, f64)> = edges
            .iter()
            .map(|&(i, j)| {
                let scale = sqrt_weight[i].min(sqrt_weight[j]);
                (-scale / sqrt_weight[i], scale / sqrt_weight[j])
            })
            .collect();

        let solution = solve(&edges, &coef, &q_energy);
        let in_x: Vec<f64> = solution
            .iter()
            .zip(&sqrt_weight)
            .map(|(u, sw)| u / sw)
            .collect();

        // `(a, b)` tuples convert to an `Observation` as *(covariate, response)*, not
        // (value, weight), so the weighted observations are built explicitly.
        let observations = targets
            .iter()
            .zip(&weights)
            .map(|(&y, &weight)| Observation {
                x: (),
                y,
                observed: (),
                weight,
            });
        let expected =
            tonic_regression_pre_sorted::<Increasing, _, _>(observations, &edges, &Average);
        for (&got, &want) in in_x.iter().zip(&expected) {
            worst = worst.max((got - want).abs());
        }
        assert!(worst < 1e-6, "seed {seed}: {in_x:?} vs {expected:?}");
    }
    println!("weighted-oracle sweep: worst deviation {worst:e}");
}

/// The oracle sweep again, but at the shipped default tolerance, where most solves finish
/// through the certified-polish path (predictor or corrector) rather than the residual
/// test. The certificate promises `||u - u*|| <= eps_abs`, so agreement with the
/// combinatorial oracle within a small multiple of `eps_abs` is exactly what acceptance
/// claims -- this is the test that would catch a certificate accepting a wrong answer.
#[test]
fn certified_finishes_match_the_exact_poset_oracle() {
    let settings = Settings {
        verbose: false,
        eps_abs: 1e-5,
        eps_rel: 1e-5,
        max_iter: 10_000,
    };
    let mut worst = 0.0f64;
    let mut instances = 0usize;
    for seed in 400..560u64 {
        let mut rng = StdRng::seed_from_u64(seed);
        let levels = rng.random_range(2..=4i64);
        let mut points: Vec<(i64, i64)> = (0..rng.random_range(2..=8usize))
            .map(|_| (rng.random_range(0..levels), rng.random_range(0..levels)))
            .collect();
        points.sort_unstable();
        points.dedup();
        let n = points.len();
        if n < 2 {
            continue;
        }
        let flat: Vec<f64> = points
            .iter()
            .flat_map(|&(a, b)| [a as f64, b as f64])
            .collect();
        let edges = derive_transitive_reduction(&flat, n, 2);
        let coef = vec![(-1.0, 1.0); edges.len()];
        let q: Vec<f64> = (0..n).map(|_| rng.random_range(-3.0..3.0f64)).collect();

        let m = edges.len();
        let constraints = Constraints {
            edges: &edges,
            coef: &coef,
            n,
        };
        let mut factor = KktFactor::new(n, &edges);
        let mut workspace = solver::Workspace::new(n, m);
        let mut polisher = Polisher::new(n, m);
        let outcome = solver::solve(
            &constraints,
            &q,
            &settings,
            &mut workspace,
            &mut factor,
            &mut polisher,
        );
        assert_eq!(outcome.status, Status::Solved, "seed {seed}");

        let targets: Vec<f64> = q.iter().map(|value| -value).collect();
        let expected = tonic_regression_pre_sorted::<Increasing, _, _>(
            targets.iter().copied(),
            &edges,
            &Average,
        );
        for (&got, &want) in workspace.x.iter().zip(&expected) {
            worst = worst.max((got - want).abs());
        }
        assert!(
            worst < 5.0 * settings.eps_abs,
            "seed {seed}: {:?} vs {expected:?}",
            workspace.x
        );
        instances += 1;
    }
    assert!(instances >= 120, "only {instances} usable instances");
    println!("default-tolerance sweep: {instances} instances, worst deviation {worst:e}");
}

/// The multiplier peel against hand-derived values, including the split signal.
#[test]
fn peel_multipliers_hand_cases() {
    let edges = [(0usize, 1usize)];
    let coef = [(-1.0, 1.0)];
    let a = Constraints {
        edges: &edges,
        coef: &coef,
        n: 2,
    };
    let raw = [0.0, 0.0];
    let mut p = Polisher::new(2, 1);

    // q = (-1, 1): the unconstrained answer (1, -1) violates the edge, both pool to 0,
    // and stationarity u + q = A^T lambda gives lambda = 1 (the Hildreth test's value).
    let report = polish::polish(&a, &[-1.0, 1.0], |_| true, &raw, &mut p);
    assert!(report.feasible);
    let peel = p.peel_multipliers(&a, &[-1.0, 1.0], None);
    assert_eq!(peel.negative, 0);
    assert_eq!(peel.worst_row, usize::MAX);
    assert!(
        (p.multipliers[0] - 1.0).abs() < 1e-15,
        "{:?}",
        p.multipliers
    );

    // q = (-1, -2): the unconstrained answer (1, 2) is already feasible, so pooling the
    // edge anyway prices it at lambda = -1/2 -- the dual saying "split here".
    let report = polish::polish(&a, &[-1.0, -2.0], |_| true, &raw, &mut p);
    assert!(report.feasible);
    let peel = p.peel_multipliers(&a, &[-1.0, -2.0], None);
    assert_eq!(peel.negative, 1);
    assert_eq!(peel.worst_row, 0);
    assert_eq!(p.multipliers[0], 0.0, "negatives must be clamped");

    // The same split with the row blocked: still counted, no longer nominated.
    let blocked = [true];
    let peel = p.peel_multipliers(&a, &[-1.0, -2.0], Some(&blocked));
    assert_eq!(peel.negative, 1);
    assert_eq!(peel.worst_row, usize::MAX);
}

/// The tuning constants are numerics, not style: pin them so a future edit has to be
/// deliberate and is visible in a diff.
#[test]
fn tuning_constants_are_pinned() {
    assert_eq!(solver::SIGMA, 1e-6);
    assert_eq!(solver::ALPHA, 1.6);
    assert_eq!(solver::RHO_INITIAL, 0.1);
    assert_eq!(solver::RHO_MIN, 1e-6);
    assert_eq!(solver::RHO_MAX, 1e6);
    assert_eq!(solver::ADAPTIVE_RHO_TOLERANCE, 5.0);
    assert_eq!(solver::CHECK_INTERVAL, 10);
    assert_eq!(solver::F32_ERROR_PER_COND, 7.2e-9);
    assert_eq!(solver::F32_ERROR_MARGIN, 0.1);
}
