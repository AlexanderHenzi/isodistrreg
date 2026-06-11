//! Zero-weight observations are dropped during preprocessing, before any algorithm runs.
//!
//! Contract (see `IsotonicDistributionalRegressionFit::fit`): a fit with zero-weight
//! observations equals the fit on the positive-weight subsample EXACTLY — same threshold
//! grid, same covariate set, bitwise-identical CDF values (filtering happens before any
//! arithmetic). With every weight zero, the result is the empty fit (sub-CDF ≡ 0).
//!
//! Each test plants zero-weight observations that would, if they leaked through, change
//! the output shape: a unique response value (would add a threshold) and a unique
//! covariate whose observations all have weight zero (would add a grid point).

use isodistrreg::total_order::{Config, Fit};
use isodistrreg::{Increasing, IsotonicDistributionalRegressionFit, NoProgress, StochasticOrder};

fn fit_total(
    x: &[f64],
    y: &[f64],
    observed: Option<&[bool]>,
    w: &[f64],
    order: StochasticOrder,
    decreasing: bool,
) -> Fit<f64, f64> {
    Fit::<f64, f64>::fit(
        x,
        y,
        observed,
        Some(w),
        Increasing,
        order,
        decreasing,
        Config,
        &NoProgress,
    )
    .unwrap()
}

/// The base data set, plus the same data with interleaved zero-weight observations:
/// one at a fresh response value (y = 1.5), one at a fresh covariate (x = 2.5), and one
/// duplicating an existing point. `observed` patterns are chosen by the caller.
struct Padded {
    x: Vec<f64>,
    y: Vec<f64>,
    w: Vec<f64>,
    observed: Vec<bool>,
}

fn padded(observed_base: &[bool]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<bool>, Padded) {
    let x_base = vec![1.0, 1.0, 2.0, 3.0, 3.0];
    let y_base = vec![0.5, 2.0, 1.0, 0.5, 3.0];
    let w_base = vec![1.0, 0.5, 2.0, 1.0, 0.75];
    assert_eq!(observed_base.len(), x_base.len());

    // Insert zero-weight observations at positions 1, 4, 7 (fresh y, fresh x, duplicate).
    let inserts: [(f64, f64, bool); 3] = [(1.0, 1.5, true), (2.5, 0.25, false), (2.0, 1.0, true)];
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut w = Vec::new();
    let mut observed = Vec::new();
    let mut base_iter = 0;
    for i in 0..x_base.len() + inserts.len() {
        if i % 3 == 1 && i / 3 < inserts.len() {
            let (xi, yi, oi) = inserts[i / 3];
            x.push(xi);
            y.push(yi);
            w.push(0.0);
            observed.push(oi);
        } else {
            x.push(x_base[base_iter]);
            y.push(y_base[base_iter]);
            w.push(w_base[base_iter]);
            observed.push(observed_base[base_iter]);
            base_iter += 1;
        }
    }
    assert_eq!(base_iter, x_base.len());

    (
        x_base,
        y_base,
        w_base,
        observed_base.to_vec(),
        Padded { x, y, w, observed },
    )
}

fn assert_fits_identical(base: &Fit<f64, f64>, padded: &Fit<f64, f64>, label: &str) {
    assert_eq!(
        base.thresholds, padded.thresholds,
        "{label}: zero-weight observations changed the threshold grid"
    );
    assert_eq!(
        base.covariates, padded.covariates,
        "{label}: zero-weight observations changed the covariate grid"
    );
    assert_eq!(
        base.cdfs, padded.cdfs,
        "{label}: zero-weight observations changed the fitted CDF values"
    );
}

#[test]
fn zero_weight_inert_total_order_all_pipelines() {
    let observed_patterns: [(&str, Option<[bool; 5]>); 2] = [
        ("uncensored", None),
        ("censored", Some([true, false, true, true, false])),
    ];
    for (pattern_label, observed_base) in observed_patterns {
        for (order_label, order) in [
            ("sd", StochasticOrder::StochasticDominance),
            ("hro", StochasticOrder::HazardRateOrder),
        ] {
            for decreasing in [false, true] {
                let base_flags = observed_base.unwrap_or([true; 5]);
                let (xb, yb, wb, ob, pad) = padded(&base_flags);

                let base = fit_total(
                    &xb,
                    &yb,
                    observed_base.as_ref().map(|_| ob.as_slice()),
                    &wb,
                    order,
                    decreasing,
                );
                let with_zeros = fit_total(
                    &pad.x,
                    &pad.y,
                    observed_base.as_ref().map(|_| pad.observed.as_slice()),
                    &pad.w,
                    order,
                    decreasing,
                );
                assert_fits_identical(
                    &base,
                    &with_zeros,
                    &format!("{pattern_label}/{order_label}/decreasing={decreasing}"),
                );
            }
        }
    }
}

#[test]
fn all_zero_weights_yield_empty_fit() {
    for (order_label, order) in [
        ("sd", StochasticOrder::StochasticDominance),
        ("hro", StochasticOrder::HazardRateOrder),
    ] {
        for observed in [None, Some([true, false, true].as_slice())] {
            let fit = fit_total(
                &[1.0, 2.0, 3.0],
                &[1.0, 2.0, 3.0],
                observed,
                &[0.0; 3],
                order,
                false,
            );
            assert!(
                fit.thresholds.is_empty() && fit.covariates.is_empty() && fit.cdfs.is_empty(),
                "{order_label}/observed={observed:?}: all-zero weights must produce the empty \
                 fit, got thresholds {:?}",
                fit.thresholds,
            );
            // The empty fit predicts the all-zero sub-CDF.
            assert_eq!(fit.cdf(2.0).count(), 0);
        }
    }
}

#[test]
fn zero_weight_inert_plain_survival_idr() {
    // Unique x and y as the function requires; the zero-weight observation sits at a
    // fresh covariate AND fresh response, so leaking it would change the output shape.
    let run = |x: &[f64], y: &[f64], o: &[bool], w: &[f64]| {
        isodistrreg::total_order::stochastic_dominance::censored_nonrecursive(x, y, o, w).unwrap()
    };
    let base = run(
        &[1.0, 2.0, 3.0],
        &[2.0, 1.0, 3.0],
        &[true, false, true],
        &[1.0, 2.0, 1.5],
    );
    let with_zero = run(
        &[1.0, 1.5, 2.0, 3.0],
        &[2.0, 0.5, 1.0, 3.0],
        &[true, true, false, true],
        &[1.0, 0.0, 2.0, 1.5],
    );
    assert_eq!(
        base.0, with_zero.0,
        "zero-weight observation changed the output dimensions"
    );
    assert_eq!(
        base.1, with_zero.1,
        "zero-weight observation changed the fitted values"
    );
}

#[cfg(feature = "partial-order")]
mod partial_order {
    use isodistrreg::partial_order::{CovariateGroups, Fit};
    use isodistrreg::{IsotonicDistributionalRegressionFit, NoProgress, StochasticOrder};

    fn fit_po(x: &[f64], y: &[f64], observed: Option<&[bool]>, w: &[f64]) -> Fit<f64, f64> {
        Fit::<f64, f64>::fit(
            x,
            y,
            observed,
            Some(w),
            CovariateGroups::empty(2),
            StochasticOrder::StochasticDominance,
            false,
            Default::default(),
            &NoProgress,
        )
        .unwrap()
    }

    /// 2-d covariates; the zero-weight observation introduces a fresh covariate row and
    /// a fresh response. Compare predictions on the training points (QP solutions are
    /// only ~1e-5 reproducible, but identical inputs give identical solver runs — the
    /// filtered problem must be bit-identical to the base problem).
    #[test]
    fn zero_weight_inert_partial_order() {
        let x_base = [0.0, 0.0, 1.0, 1.0, 2.0, 2.0];
        let y_base = [1.0, 2.0, 1.5];
        let w_base = [1.0, 2.0, 1.0];
        // Insert a zero-weight observation with fresh covariate row and response.
        let x_pad = [0.0, 0.0, 5.0, 5.0, 1.0, 1.0, 2.0, 2.0];
        let y_pad = [1.0, 9.0, 2.0, 1.5];
        let w_pad = [1.0, 0.0, 2.0, 1.0];

        for observed in [None, Some(())] {
            let (ob, op): (Option<Vec<bool>>, Option<Vec<bool>>) = match observed {
                None => (None, None),
                Some(()) => (
                    Some(vec![true, false, true]),
                    Some(vec![true, true, false, true]),
                ),
            };
            let base = fit_po(&x_base, &y_base, ob.as_deref(), &w_base);
            let padded = fit_po(&x_pad, &y_pad, op.as_deref(), &w_pad);

            assert_eq!(base.thresholds(), padded.thresholds());
            let q_base: Vec<f32> = base.cdf(&[1.0, 1.0]).collect();
            let q_padded: Vec<f32> = padded.cdf(&[1.0, 1.0]).collect();
            assert_eq!(
                q_base, q_padded,
                "observed={observed:?}: zero-weight observation changed the fit"
            );
        }
    }

    #[test]
    fn all_zero_weights_yield_empty_fit_partial_order() {
        for observed in [None, Some([true, false].as_slice())] {
            let fit = fit_po(&[0.0, 0.0, 1.0, 1.0], &[1.0, 2.0], observed, &[0.0, 0.0]);
            assert!(
                fit.thresholds().is_empty(),
                "observed={observed:?}: all-zero weights must produce the empty fit"
            );
            assert_eq!(fit.cdf(&[0.5, 0.5]).count(), 0);
        }
    }
}
