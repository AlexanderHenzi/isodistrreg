//! `mean()` distinguishes a proper CDF from a sub-CDF by the EXACT test
//! `last CDF value < 1.0` (src/prediction.rs). That is only sound because every
//! producer pins values that are mathematically exactly 1, instead of leaving them a
//! few ulps short through f32/f64 running-weight round-off:
//!
//! - the uncensored kernels write literal 1.0 rows for the last threshold,
//! - `routines::empirical_cdf` pins its final value (total/total = 1),
//! - the censored Kaplan-Meier producers pin survival to exactly 0 whenever the last
//!   positive-weight observation (response order, events before censorings at ties) is
//!   an EVENT — the purely combinatorial condition under which some K-M factor is
//!   exactly `1 - w/w`.
//!
//! The weights below (0.1, 0.3, 0.7, ...) are deliberately not exactly representable
//! in binary so that naive running sums/subtractions drift off the exact zero.

use isodistrreg::routines::kaplan_meier;
use isodistrreg::total_order::{Config, Fit};
use isodistrreg::{
    Increasing, IsotonicDistributionalRegressionFit, NoProgress, Observation, StochasticOrder,
};

fn obs(y: f64, observed: bool, weight: f32) -> Observation<(), f64, bool, f32> {
    Observation {
        x: (),
        y,
        observed,
        weight,
    }
}

fn fit(
    x: &[f64],
    y: &[f64],
    observed: &[bool],
    w: &[f64],
    order: StochasticOrder,
) -> Fit<f64, f64> {
    Fit::<f64, f64>::fit(
        x,
        y,
        Some(observed),
        Some(w),
        Increasing,
        order,
        false,
        Config,
        &NoProgress,
    )
    .unwrap()
}

#[test]
fn kaplan_meier_completing_data_ends_exactly_at_one() {
    // Last positive-weight observation is an event => survival exactly 0.
    let observations = vec![
        obs(1.0, true, 0.1),
        obs(2.0, false, 0.3),
        obs(3.0, true, 0.45),
        obs(4.0, true, 0.7),
    ];
    let total: f32 = 0.1 + 0.3 + 0.45 + 0.7;
    let result = kaplan_meier(observations.into_iter(), total);
    assert_eq!(
        *result.last().unwrap(),
        1.0,
        "completing data must end at exactly 1.0, got {result:?}"
    );
}

#[test]
fn kaplan_meier_trailing_censoring_stays_below_one() {
    // Last positive-weight observation censored => mass remains beyond the grid.
    let observations = vec![obs(1.0, true, 0.3), obs(2.0, false, 0.1)];
    let result = kaplan_meier(observations.into_iter(), 0.4);
    assert!(
        *result.last().unwrap() < 1.0,
        "trailing censoring must keep the sub-CDF below 1, got {result:?}"
    );
}

/// SD-censored fast path: x=1 has an event then a censoring, x=2 ends in an event.
/// Cells {2} and {1,2} end in events (survival exactly 0); the isotonic min over
/// cells pulls covariate 1 to the pooled cell, so BOTH covariates have exact unit
/// mass at the last threshold and finite means.
#[test]
fn sd_censored_completing_fit_has_exact_unit_mass_and_finite_mean() {
    let f = fit(
        &[1.0, 1.0, 2.0],
        &[2.0, 1.0, 3.0],
        &[false, true, true],
        &[0.1, 0.3, 0.7],
        StochasticOrder::StochasticDominance,
    );
    assert_eq!(f.thresholds, vec![1.0, 3.0]);
    let n_threshold = f.thresholds.len();
    for (c, &x) in f.covariates.clone().iter().enumerate() {
        let last = f.cdfs[c * n_threshold + n_threshold - 1];
        assert_eq!(last, 1.0, "covariate {x}: last CDF must be exactly 1.0");
        let mean = f.mean(x);
        assert!(
            mean.is_finite(),
            "covariate {x}: mean must be finite, got {mean}"
        );
    }
}

/// Same data through the hazard-rate censored kernel: covariate 1's at-risk mass is
/// pinned to exactly 0 when its last (censored) observation leaves, covariate 2's raw
/// estimator to exactly 0 at its final event, so the pooled ratio is exactly 0 and
/// both survivals end at exactly 0.
#[test]
fn hazard_censored_completing_fit_has_exact_unit_mass_and_finite_mean() {
    let f = fit(
        &[1.0, 1.0, 2.0],
        &[2.0, 1.0, 3.0],
        &[false, true, true],
        &[0.1, 0.3, 0.7],
        StochasticOrder::HazardRateOrder,
    );
    assert_eq!(f.thresholds, vec![1.0, 3.0]);
    let n_threshold = f.thresholds.len();
    for (c, &x) in f.covariates.clone().iter().enumerate() {
        let last = f.cdfs[c * n_threshold + n_threshold - 1];
        assert_eq!(last, 1.0, "covariate {x}: last CDF must be exactly 1.0");
        let mean = f.mean(x);
        assert!(
            mean.is_finite(),
            "covariate {x}: mean must be finite, got {mean}"
        );
    }
}

/// Single-covariate hazard path (the Kaplan-Meier curve filtered to event
/// thresholds): the final event completes the distribution.
#[test]
fn hazard_censored_single_covariate_completing_mean_is_finite() {
    let f = fit(
        &[5.0, 5.0, 5.0],
        &[1.0, 2.0, 3.0],
        &[true, false, true],
        &[0.1, 0.3, 0.7],
        StochasticOrder::HazardRateOrder,
    );
    assert_eq!(f.thresholds, vec![1.0, 3.0]);
    assert_eq!(*f.cdfs.last().unwrap(), 1.0);
    assert!(f.mean(5.0).is_finite());
}

/// The exact gate must keep flagging genuine sub-CDFs: trailing censored mass means
/// the mean is undefined (NaN), for both censored kernels.
#[test]
fn censored_sub_cdf_mean_stays_nan() {
    for (name, order) in [
        ("sd", StochasticOrder::StochasticDominance),
        ("hazard", StochasticOrder::HazardRateOrder),
    ] {
        let f = fit(&[1.0, 1.0], &[1.0, 2.0], &[true, false], &[0.3, 0.1], order);
        let mean = f.mean(1.0);
        assert!(
            mean.is_nan(),
            "{name}: trailing censored mass must give a NaN mean, got {mean}"
        );
    }
}
