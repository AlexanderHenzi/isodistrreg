use crossbeam_channel::RecvTimeoutError;
use extendr_api::scalar::Rbool;
use extendr_api::symbol::class_symbol;
use extendr_api::{
    Attributes, IntoRobj, List, Nullable, Operators, RColumn, RMatrix, Rinternals, Robj, extendr,
    extendr_module,
};
use isodistrreg::functionals::{ClippingWrapper, KaplanMeier};
use isodistrreg::partial_order::{CovariateGroups, PredictionWorkspace};
use isodistrreg::{
    CovariateInterpolator, Decreasing, Increasing, NoProgress, Parallel, ProgressTracker, Serial,
    StochasticOrder, quantile, subagging,
};
use isodistrreg::{IsotonicDistributionalRegressionFit, partial_order, total_order};
use itertools::{Either, izip};
use std::iter::repeat;
use std::str::FromStr;

mod logging;

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod isodistrreg;
    // core methods
    impl IDR;
    // additional methods
    fn survival_isotonic_distributional_regression_threshold;
    fn plain_survival_isotonic_distributional_regression;
    fn plain_survival_isotonic_distributional_regression_threshold;
    // utility
    fn isotonic_regression_impl;
}

#[allow(clippy::upper_case_acronyms)]
#[extendr]
enum IDR {
    Partial(subagging::Fit<partial_order::Fit<f64, f64>>),
    Total(subagging::Fit<total_order::Fit<f64, f64>>),
}

#[allow(non_snake_case)]
#[extendr]
impl IDR {
    /// Wrapper around the Rust function. No careful input parsing here, we expect that to have
    /// happened on the R side.
    ///
    /// @description Internal method that is used by R code to pass into Rust code.
    /// @param y Outcome / response values (might be observed or censored)
    /// @param X Flattened array of covariate values.
    /// @param y_observed Nullable array of the same length as the responses indicating censoring
    ///   (TRUE represents observed, FALSE represents right-censored), default is all observed
    /// @param sample_weight Nullable array of the same length as the responses with nonnegative weights.
    /// @param x_order List of covariate group descriptions.
    /// @param y_order String representing the kind of stochastic order of the response.
    /// @param decreasing Whether the responses are decreasing with the covariate (not is increasing).
    /// @param settings List of solver settings for partially ordered covariates.
    /// @param seed Integer seed for the random number generator. Only relevant when (su)bagging
    ///   is active (i.e. when `subsamples` is set).
    /// @param n_jobs Number of worker threads used to fit the individual subsamples in parallel.
    ///   Only relevant when (su)bagging is active. Default is 1 (serial execution).
    /// @param show_progress Whether to display a progress bar while fitting. Default is false.
    #[allow(clippy::too_many_arguments)]
    fn fit(
        y: &[f64],
        X: &[f64],
        y_observed: Nullable<&[Rbool]>,
        sample_weight: Nullable<&[f64]>,
        x_order: List,
        y_order: &str,
        decreasing: bool,
        subsamples: Nullable<usize>,
        subsample_size: Nullable<f64>,
        replace: Nullable<bool>,
        settings: Nullable<List>,
        seed: Nullable<f64>,
        n_jobs: Nullable<usize>,
        show_progress: Nullable<bool>,
    ) -> Self {
        assert!(!y.is_empty(), "'y' must not be empty");
        assert_eq!(
            X.len() % y.len(),
            0,
            "Covariates and responses shapes are not compatible",
        );
        let x_dimension = X.len() / y.len();

        let maybe_observed = y_observed
            .map(|o| o.iter().map(Rbool::to_bool).collect::<Vec<_>>())
            .into_option();
        let maybe_sample_weight = match sample_weight {
            Nullable::NotNull(slice) => Some(slice),
            Nullable::Null => None,
        };

        // may be empty
        let x_order_parsed = CovariateGroups::parse(
            x_order.into_iter().map(|(name, value)| {
                let members = value
                    .as_integer_vector()
                    .expect("covariate group indices must be integers")
                    .into_iter()
                    .map(|i| {
                        assert!(
                            i != i32::MIN && i >= 1,
                            "covariate group indices must be positive integers without NA",
                        );
                        // R indexing starts with 1
                        i as usize - 1
                    });
                (name, members)
            }),
            x_dimension,
        )
        .unwrap_or_else(|e| panic!("{e}"));

        let y_order_parsed = StochasticOrder::from_str(y_order).unwrap_or_else(|e| panic!("{e}"));
        let config = subagging::Config::parse(
            subsamples.into_option(),
            subsample_size.into_option().map(Either::Right),
            replace.into_option().unwrap_or(false),
            y.len(),
            seed.into_option().map(|s| s as u64),
            n_jobs.into_option().unwrap_or(1),
        )
        .unwrap();

        let mut osqp_settings = osqp::Settings::default();
        if let Some(user_settings) = settings.into_option() {
            osqp_settings = parse_settings(user_settings, osqp_settings);
        }

        // The fit runs on a worker thread when a progress bar is enabled, so the main thread is
        // free to drain progress messages and render through `RTerm` (which calls into R).
        // Worker threads must never touch R themselves.
        let progress_pair = match show_progress {
            Nullable::NotNull(true) => Some(logging::ProgressPump::new()),
            _ => None,
        };
        let progress: &dyn ProgressTracker = match &progress_pair {
            Some((tracker, _)) => tracker,
            None => &NoProgress,
        };

        let run_fit = move || -> IDR {
            match x_dimension {
                0 => panic!("Covariate of empty dimension"),
                1 => {
                    let fit = subagging::Fit::<total_order::Fit<f64, f64>>::fit(
                        X,
                        y,
                        maybe_observed.as_deref(),
                        maybe_sample_weight,
                        Default::default(),
                        y_order_parsed,
                        decreasing,
                        (config, Default::default()),
                        progress,
                    )
                    .unwrap_or_else(|e| panic!("invalid input: {e}"));
                    IDR::Total(fit)
                }
                _ => {
                    let fit = subagging::Fit::<partial_order::Fit<f64, f64>>::fit(
                        X,
                        y,
                        maybe_observed.as_deref(),
                        maybe_sample_weight,
                        x_order_parsed,
                        y_order_parsed,
                        decreasing,
                        (config, partial_order::Config { osqp_settings }),
                        progress,
                    )
                    .unwrap_or_else(|e| panic!("invalid input: {e}"));
                    IDR::Partial(fit)
                }
            }
        };

        let result = match progress_pair.as_ref() {
            Some((_, pump)) => std::thread::scope(|s| {
                let handle = s.spawn(run_fit);
                loop {
                    match pump.recv() {
                        Ok(msg) => {
                            pump.apply(msg);
                            pump.drain();
                        }
                        Err(RecvTimeoutError::Timeout) => {
                            if handle.is_finished() {
                                break;
                            }
                        }
                        Err(RecvTimeoutError::Disconnected) => break,
                    }
                }
                match handle.join() {
                    Ok(idr) => idr,
                    Err(payload) => std::panic::resume_unwind(payload),
                }
            }),
            None => run_fit(),
        };

        if let Some((_, pump)) = &progress_pair {
            pump.finish();
        }

        result
    }

    /// Predict a conditional CDF for the covariates provided in `data`.
    ///
    /// @param data Covariate array (flattened) for which to predict CDFs, (covariate-major) layout
    ///   is inferred using this fit's covariate dimension.
    fn cdf(&self, X: &[f64]) -> RMatrix<f64> {
        // R's REALSXP is always 64-bit double — there is no float32 in R. We promote element-wise
        // at the binding boundary. This copy was already mandatory because of how
        // `RMatrix::new_matrix` works, so there's no extra cost.
        match self {
            IDR::Partial(fit) => {
                let flat: Vec<_> = X
                    .chunks_exact(fit.covariate_groups.dimension)
                    .flat_map(|c| fit.cdf(c))
                    .collect();
                let n = X.len() / fit.covariate_groups.dimension;
                RMatrix::new_matrix(n, fit.thresholds.len(), |r, c| {
                    flat[r * fit.thresholds.len() + c] as f64
                })
            }
            IDR::Total(fit) => {
                let flat: Vec<_> = X.iter().flat_map(|&c| fit.cdf(c)).collect();
                RMatrix::new_matrix(X.len(), fit.thresholds.len(), |r, c| {
                    flat[r * fit.thresholds.len() + c] as f64
                })
            }
        }
    }

    fn cdf_at(&self, X: &[f64], y: &[f64]) -> RMatrix<f64> {
        match self {
            IDR::Partial(fit) => {
                let mut workspace = PredictionWorkspace::new();
                let interpolations: Vec<_> = X
                    .chunks_exact(fit.covariate_groups.dimension)
                    .map(|c| fit.interpolate_covariate_with_workspace(c, &mut workspace))
                    .collect();
                let coordinates: Vec<_> =
                    y.iter().map(|&t| fit.get_response_coordinate(t)).collect();

                let n_covariate = interpolations.len();
                RMatrix::new_matrix(n_covariate, y.len(), |r, c| {
                    // f32 → f64 promotion at the R boundary; see `cdf`.
                    interpolations[r].interpolate(coordinates[c]) as f64
                })
            }
            IDR::Total(fit) => {
                let interpolations: Vec<_> =
                    X.iter().map(|&c| fit.interpolate_covariate(c)).collect();
                let coordinates: Vec<_> =
                    y.iter().map(|&t| fit.get_response_coordinate(t)).collect();
                let n_covariate = interpolations.len();
                RMatrix::new_matrix(n_covariate, y.len(), |r, c| {
                    interpolations[r].interpolate(coordinates[c]) as f64
                })
            }
        }
    }

    fn quantile(&self, X: &[f64], probability: &[f64]) -> RMatrix<f64> {
        match self {
            IDR::Partial(fit) => {
                let mut workspace = PredictionWorkspace::new();
                let interpolations: Vec<_> = X
                    .chunks_exact(fit.covariate_groups.dimension)
                    .map(|c| fit.interpolate_covariate_with_workspace(c, &mut workspace))
                    .collect();
                let n_covariate = interpolations.len();
                RMatrix::new_matrix(n_covariate, probability.len(), |r, c| {
                    quantile(
                        &interpolations[r],
                        probability[c] as f32,
                        false,
                        &fit.thresholds,
                    )
                })
            }
            IDR::Total(fit) => {
                let interpolations: Vec<_> =
                    X.iter().map(|&c| fit.interpolate_covariate(c)).collect();
                let n_covariate = interpolations.len();
                RMatrix::new_matrix(n_covariate, probability.len(), |r, c| {
                    quantile(
                        &interpolations[r],
                        probability[c] as f32,
                        false,
                        &fit.thresholds,
                    )
                })
            }
        }
    }

    fn dimension(&self) -> usize {
        match self {
            IDR::Partial(fit) => fit.covariate_groups.dimension,
            IDR::Total(_) => 1,
        }
    }

    fn thresholds(&self) -> &[f64] {
        match self {
            IDR::Partial(fit) => &fit.thresholds,
            IDR::Total(fit) => &fit.thresholds,
        }
    }

    fn diagnostic(&self) -> List {
        match self {
            IDR::Partial(subagging::Fit { fits, .. }) => {
                let worst = fits
                    .iter()
                    .map(
                        |&partial_order::Fit {
                             quality_indicators, ..
                         }| quality_indicators,
                    )
                    .reduce(|a, b| partial_order::QualityIndicators {
                        precision: a.precision.max(b.precision),
                        convergence_fraction: a.convergence_fraction.min(b.convergence_fraction),
                    })
                    .unwrap();
                List::from_pairs([
                    ("precision", worst.precision.into_robj()),
                    (
                        "convergence_fraction",
                        worst.convergence_fraction.into_robj(),
                    ),
                ])
            }
            IDR::Total(subagging::Fit { fits, .. }) => {
                let worst = fits
                    .iter()
                    .map(
                        |&total_order::Fit {
                             quality_indicators, ..
                         }| quality_indicators,
                    )
                    .reduce(|a, b| total_order::QualityIndicators {
                        epsilon: a.epsilon.max(b.epsilon),
                    })
                    .unwrap();
                List::from_pairs([("epsilon", worst.epsilon.into_robj())])
            }
        }
    }
}

fn parse_settings(user_settings: List, mut existing: osqp::Settings) -> osqp::Settings {
    // Missing keys (`$` yields NULL) keep the defaults; present-but-invalid values must
    // produce a clear error here rather than an unwrap panic or an OSQP setup failure.
    if let Ok(value) = user_settings.dollar("verbose")
        && !value.is_null()
    {
        let v = value
            .as_bool()
            .expect("'pars$verbose' must be TRUE or FALSE");
        existing = existing.verbose(v);
    }
    if let Ok(value) = user_settings.dollar("eps_abs")
        && !value.is_null()
    {
        let v = value
            .as_real()
            .filter(|v| v.is_finite() && *v >= 0.0)
            .expect("'pars$eps_abs' must be a finite non-negative number");
        existing = existing.eps_abs(v);
    }
    if let Ok(value) = user_settings.dollar("eps_rel")
        && !value.is_null()
    {
        let v = value
            .as_real()
            .filter(|v| v.is_finite() && *v >= 0.0)
            .expect("'pars$eps_rel' must be a finite non-negative number");
        existing = existing.eps_rel(v);
    }
    if let Ok(value) = user_settings.dollar("max_iter")
        && !value.is_null()
    {
        let v = value
            .as_integer()
            .map(f64::from)
            .or_else(|| value.as_real())
            .filter(|v| v.fract() == 0.0 && *v >= 1.0 && *v <= f64::from(u32::MAX))
            .expect("'pars$max_iter' must be a positive integer");
        existing = existing.max_iter(v as u32);
    }

    existing
}

/// Compute the isotonic regression for the mean for totally ordered covariates.
///
/// Internal wrapper; user-facing input validation happens in the R function
/// `isotonic_regression()` (R/modeling.R). The asserts here are backstops with clear
/// messages for anyone calling the wrapper directly.
///
/// @description Internal method that is used by R code to pass into Rust code.
/// @param y Double vector of response values.
/// @param X Double vector of covariate values, or NULL if responses are pre-sorted.
/// @param weights Double vector of non-negative weights, or NULL for equal weights.
/// @param decreasing Bool indicating direction (default FALSE is increasing, TRUE is decreasing).
/// @returns Numeric vector of isotonic fitted means.
#[allow(non_snake_case)]
#[extendr]
fn isotonic_regression_impl(
    y: &[f64],
    #[extendr(default = "NULL")] X: Nullable<&[f64]>,
    #[extendr(default = "NULL")] weights: Nullable<&[f64]>,
    #[extendr(default = "FALSE")] decreasing: bool,
) -> Vec<f64> {
    let n = y.len();
    assert!(
        y.iter().all(|v| v.is_finite()),
        "'y' must contain only finite values",
    );

    macro_rules! tonic {
        ($func:ident, $data:expr) => {
            match decreasing {
                false => total_order::$func::<Increasing, f64, _>($data),
                true => total_order::$func::<Decreasing, f64, _>($data),
            }
        };
    }

    if let Nullable::NotNull(c) = X {
        assert_eq!(c.len(), n, "'X' and 'y' must have equal length");
        assert!(
            c.iter().all(|v| v.is_finite()),
            "'X' must contain only finite values",
        );
    }
    if let Nullable::NotNull(w) = weights {
        assert_eq!(w.len(), n, "'weights' and 'y' must have equal length");
        assert!(
            w.iter().all(|v| v.is_finite() && *v >= 0.0),
            "'weights' must contain only finite non-negative values",
        );
        assert!(
            n == 0 || w.iter().any(|&v| v > 0.0),
            "at least one weight must be positive",
        );
    }

    let partition_iterator = match (X, weights) {
        (Nullable::NotNull(c), Nullable::NotNull(w)) => {
            let data = izip!(c.iter().copied(), y.iter().copied(), w.iter().copied());
            tonic!(tonic_regression, data)
        }
        (Nullable::NotNull(c), Nullable::Null) => {
            let data = izip!(c.iter().copied(), y.iter().copied(), repeat(1.0));
            tonic!(tonic_regression, data)
        }
        (Nullable::Null, Nullable::NotNull(w)) => {
            let data = y.iter().copied().zip(w.iter().copied());
            tonic!(tonic_regression_pre_sorted, data)
        }
        (Nullable::Null, Nullable::Null) => {
            let data = y.iter().copied().zip(repeat(1.0));
            tonic!(tonic_regression_pre_sorted, data)
        }
    };

    partition_iterator.collect()
}

/// Compute an isotonic regression of a probability given by Kaplan-Meier estimators (= one
/// threshold of S-IDR).
///
/// This method is provided for completeness and to make the S-IDR publication easier to reproduce.
///
/// @param threshold Double of the response value at which to compute the IDR solution.
/// @param X Double vector of totally ordered covariates.
/// @param y Double vector of response values.
/// @param y_observed Integer vector: 1 for observed, 0 for censored.
/// @param weights Double vector of non-negative weights. All vectors must have equal length.
/// @param decreasing Bool indicating direction (decreasing is a CIDR threshold).
/// @returns Numeric vector. The i'th entry gives the fitted CDF at covariate i.
/// @examples
/// survival_isotonic_distributional_regression_threshold(3.5, as.double(1:4), c(2, 1, 4, 3),
///   as.integer(c(1, 0, 0, 0)), rep(1.0, 4))
/// @export
#[allow(non_snake_case)]
#[extendr]
fn survival_isotonic_distributional_regression_threshold(
    threshold: f64,
    X: &[f64],
    y: &[f64],
    y_observed: &[i32],
    weights: &[f64],
    #[extendr(default = "FALSE")] decreasing: bool,
) -> RColumn<f64> {
    assert!(!y.is_empty(), "'y' must not be empty");
    assert_eq!(X.len(), y.len(), "'X' and 'y' must have equal length");
    assert_eq!(
        y_observed.len(),
        y.len(),
        "'y_observed' and 'y' must have equal length",
    );
    assert_eq!(
        weights.len(),
        y.len(),
        "'weights' and 'y' must have equal length",
    );
    assert!(threshold.is_finite(), "'threshold' must be finite");
    assert!(
        X.iter().chain(y).chain(weights).all(|v| v.is_finite()),
        "'X', 'y' and 'weights' must contain only finite values",
    );
    assert!(
        weights.iter().all(|&w| w >= 0.0),
        "'weights' must be non-negative",
    );
    assert!(
        y_observed.iter().all(|&c| c != i32::MIN),
        "'y_observed' must not contain NA",
    );

    let data = izip!(
        X.iter().copied(),
        y.iter().copied(),
        y_observed.iter().map(|&c| c > 0),
        weights.iter().copied(),
    );
    let functional = ClippingWrapper::new(KaplanMeier::new(threshold));
    let result = if decreasing {
        total_order::functionals::algorithm::<Decreasing, _, _>(data, &functional)
    } else {
        total_order::functionals::algorithm::<Increasing, _, _>(data, &functional)
    };

    // KM/PAVA now returns f32 (matching the rest of the censored stack); widen to f64 for
    // the R column constructor.
    RColumn::new_column(result.len(), |i| f64::from(result[i]))
}

/// Compute the plain survival IDR for totally ordered co-variates under the hazard rate order
/// assumption for censored data. Note that the hazard rate order is needed for the method to be
/// consistent, but that the returned solution likely doesn't satisfy this assumption.
///
/// This method is provided for completeness and to make the S-IDR publication easier to reproduce,
/// but is not recommended for use by practitioners due to prohibitive computational costs and
/// potential inconsistency in case the distributions are not hazard rate ordered.
///
/// @description
/// Computes the plain survival isotonic distributional regression (plain survival IDR) under the
/// hazard rate order assumption for totally ordered covariates when some responses are right
/// censored. Returns the fitted cumulative distribution evaluated at each response threshold for
/// each unique covariate.
///
/// @param X Double vector of scalar covariates.
/// @param y Double vector of response values.
/// @param y_observed Integer vector: 1 for observed, 0 for censored.
/// @param weights Double vector of non-negative weights. All vectors must have equal length.
/// @returns Numeric matrix. The (i, j) entry gives the fitted CDF at (unique) response j for
///   (unique) covariate i.
/// @examples
/// plain_survival_isotonic_distributional_regression(as.double(1:4), c(2, 1, 4, 3),
///   as.integer(c(0, 1, 0, 1)), c(1, 2, 1, 1))
/// @export
#[allow(non_snake_case)]
#[extendr]
fn plain_survival_isotonic_distributional_regression(
    X: &[f64],
    y: &[f64],
    y_observed: &[i32],
    weights: &[f64],
) -> RMatrix<f64> {
    let observed: Vec<_> = y_observed.iter().map(|&c| c > 0).collect();
    let ((n_covariate, n_threshold), result) =
        total_order::stochastic_dominance::censored_nonrecursive(X, y, &observed, weights)
            .unwrap_or_else(|e| panic!("invalid input: {e}"));

    RMatrix::new_matrix(n_covariate, n_threshold, |covariate, response| {
        result[response * n_covariate + covariate]
    })
}

/// Compute an isotonic regression of a probability given by Kaplan-Meier estimators (= one
/// threshold of the plain survival IDR data). Note that this is the non-recursive
/// version that is consistent only under the hazard rate order assumption.
///
/// @description
/// Computes a single threshold of the plain survival isotonic distributional regression (plain
/// survival IDR) under the hazard rate order assumption for totally ordered covariates when some
/// responses are right-censored. Returns the fitted cumulative distribution evaluated at the
/// threshold for each unique covariate.
///
/// @param threshold Double of the response value at which to compute the IDR solution.
/// @param X Double vector of totally ordered covariates.
/// @param y Double vector of response values.
/// @param y_observed Integer vector: 1 for observed, 0 for censored.
/// @param weights Double vector of non-negative weights. All vectors must have equal length.
/// @param decreasing Bool indicating direction (decreasing is a CIDR threshold).
/// @param parallel Bool indicating whether to use multiple cores.
/// @returns Numeric vector. The i'th entry gives the fitted CDF at covariate i.
/// @examples
/// plain_survival_isotonic_distributional_regression_threshold(3.5, as.double(1:4), c(2, 1, 4, 3),
///   as.integer(c(1, 0, 0, 0)), rep(1.0, 4), decreasing = TRUE)
/// @export
#[allow(non_snake_case)]
#[extendr]
fn plain_survival_isotonic_distributional_regression_threshold(
    threshold: f64,
    X: &[f64],
    y: &[f64],
    y_observed: &[i32],
    weights: &[f64],
    #[extendr(default = "FALSE")] decreasing: bool,
    #[extendr(default = "FALSE")] parallel: bool,
) -> RColumn<f64> {
    let observed: Vec<_> = y_observed.iter().map(|&c| c > 0).collect();
    let result = match (decreasing, parallel) {
        (true, false) => total_order::stochastic_dominance::censored_nonrecursive_single::<
            Decreasing,
            Serial,
        >(threshold, X, y, &observed, weights),
        (true, true) => total_order::stochastic_dominance::censored_nonrecursive_single::<
            Decreasing,
            Parallel,
        >(threshold, X, y, &observed, weights),
        (false, false) => total_order::stochastic_dominance::censored_nonrecursive_single::<
            Increasing,
            Serial,
        >(threshold, X, y, &observed, weights),
        (false, true) => total_order::stochastic_dominance::censored_nonrecursive_single::<
            Increasing,
            Parallel,
        >(threshold, X, y, &observed, weights),
    }
    .unwrap_or_else(|e| panic!("invalid input: {e}"));

    RColumn::new_column(result.len(), |i| result[i])
}
