use crate::park::ParkFitResult;
use isodistrreg::functionals::Average;
use isodistrreg::partial_order::{
    CovariateGroups, OrderingInfo, PredictionWorkspace, QualityIndicators,
    derive_transitive_reduction, preprocess_uncensored,
};
use isodistrreg::routines::{argsort_unstable_by, lexicographic_cmp};
use isodistrreg::{CovariateInterpolator, Observation, ResponseCoordinate};
use isodistrreg::{Decreasing, Direction, Increasing, StochasticOrder};
use isodistrreg::{Error, quantile};
use isodistrreg::{IsotonicDistributionalRegressionFit, NoProgress, ProgressTracker, subagging};
use isodistrreg::{partial_order, total_order};
use itertools::EitherOrBoth::{Both, Left, Right};
use itertools::{Either, EitherOrBoth, izip};
use kdam::{Bar, BarExt, tqdm};
use numpy::ndarray::{
    Array, Array2, ArrayD, ArrayView, ArrayView1, ArrayView2, ArrayViewD, ArrayViewMut1, Axis,
    Dimension, Zip,
};
use numpy::{
    AllowTypeChange, Element, IntoPyArray, IxDyn, PyArray, PyArray1, PyArray2, PyArrayDescrMethods,
    PyArrayDyn, PyArrayLike, PyArrayLike1, PyArrayLike2, PyArrayLikeDyn, PyArrayMethods,
    PyUntypedArray, PyUntypedArrayMethods, dtype,
};
use pyo3::exceptions::{PyException, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple, PyType};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::mem;
use std::sync::{Mutex, OnceLock};

/// Isotonic distributional regression (IDR) model.
///
/// Fits a nonparametric distributional regression model that predicts full
/// conditional distributions (not just means) of a response variable given
/// numeric or ordinal covariates, under the assumption that the response
/// distribution increases stochastically with the covariates.
///
/// The model is fitted by calling the constructor:
///
///     model = IDR(y, X)
///
/// For univariate covariates, a total order is used. For multivariate
/// covariates, a partial (componentwise) order is used by default.
///
/// The fitted model supports prediction of conditional means, full CDFs,
/// CDF values at specific thresholds, and quantiles.
///
/// Parameters
/// ----------
/// y : array_like, shape (n,)
///     Response values (observations).
/// X : array_like, shape (n,) or (n, d)
///     Covariate values. A 1-D array is treated as a single covariate
///     with a total order. A 2-D array with d > 1 columns uses a partial
///     (componentwise) order by default.
/// y_observed : array_like of bool, shape (n,), optional
///     If given, indicates right-censored observations (True = observed,
///     False = censored). If not provided, interpreted as all observed.
/// sample_weight : array_like, shape (n,), optional
///     Non-negative observation weights.
/// X_order : list of (str, list of int), optional
///     Specifies the partial order on the covariate space. Each entry is
///     a tuple ``(kind, column_indices)`` where *kind* is one of
///     ``"comp"`` (componentwise), ``"sd"`` (stochastic dominance), or
///     ``"icx"`` (increasing convex). Only used for multivariate
///     covariates.
/// y_order : str, optional
///     Stochastic order on the response: ``"sd"`` for stochastic dominance
///     (default) or ``"hazard"`` for hazard rate order.
/// decreasing : bool, optional
///     If True, fit a decreasing (antitone) model instead of an
///     increasing (isotone) one. Default is False.
/// subsamples : int, optional
///     Number of (su)bagging iterations. When set, the model is fitted on
///     random subsamples and the results are averaged.
/// subsample_size : int or float, optional
///     Size of each subsample. An int is an absolute count; a float in
///     (0, 1] is a fraction of the training set.
/// replace : bool, optional
///     If True, sample with replacement (bootstrap). Default is False.
/// settings : dict, optional
///     Solver settings. For multivariate covariates this may contain an
///     ``"osqp_settings"`` dict with keys ``"verbose"`` (bool),
///     ``"eps_abs"`` (float), ``"eps_rel"`` (float), ``"max_iter"``
///     (int). For univariate covariates this may contain ``"epsilon"``
///     (float), cdf values smaller than this may be rounded.
/// seed : int, optional
///     Random seed for reproducibility. Only relevant when (su)bagging is
///     active (i.e. when ``subsamples`` is set).
/// n_jobs : int, optional
///     Number of worker threads used to fit the individual subsamples in
///     parallel. Only relevant when (su)bagging is active. Default is 1
///     (serial execution).
/// progress : bool, optional
///     If True, display a progress bar while fitting. Default is False.
///     In Jupyter notebooks the bar renders as an ipywidget; in terminals
///     it renders as a tqdm-style ANSI bar on stderr.
///
/// Raises
/// ------
/// ValueError
///     If inputs have incompatible shapes or invalid values.
#[allow(clippy::upper_case_acronyms)]
#[pyclass(module = "isodistrreg._core")]
struct IDR {
    inner: Fit,
}

#[derive(Deserialize, Serialize)]
enum Fit {
    Partial(subagging::Fit<partial_order::Fit>),
    Total {
        fit: subagging::Fit<total_order::Fit>,
        squeeze: bool,
    },
}

impl IDR {
    /// Returns (order, n_covariates, covariate_dim, n_thresholds, n_subsamples, increasing).
    fn summary(&self) -> (&str, usize, usize, usize, usize, bool) {
        match &self.inner {
            Fit::Total { fit, .. } => (
                "total",
                fit.covariates.len(),
                1,
                fit.thresholds.len(),
                fit.fits.len(),
                fit.fits.first().is_none_or(|f| f.increasing),
            ),
            Fit::Partial(fit) => {
                let dim = fit.covariate_groups.dimension;
                (
                    "partial",
                    fit.covariates.len() / dim,
                    dim,
                    fit.thresholds.len(),
                    fit.fits.len(),
                    fit.fits.first().is_none_or(|f| f.increasing),
                )
            }
        }
    }
}

#[allow(non_snake_case)]
#[pymethods]
impl IDR {
    /// Fit an isotonic distributional regression model (see class documentation).
    #[allow(clippy::too_many_arguments)]
    #[new]
    #[pyo3(
        signature = (
            y,
            X,
            y_observed=None,
            sample_weight=None,
            X_order=None,
            y_order="sd",
            decreasing=false,
            subsamples=None,
            subsample_size=None,
            replace=false,
            settings=None,
            seed=None,
            n_jobs=1,
            progress=false,
        )
    )]
    fn fit(
        py: Python,
        y: PyArrayLike1<f64, AllowTypeChange>,
        X: PyArrayLikeDyn<f64, AllowTypeChange>,
        y_observed: Option<PyArrayLike1<bool, AllowTypeChange>>,
        sample_weight: Option<PyArrayLike1<f64, AllowTypeChange>>,
        X_order: Option<Vec<(String, Vec<usize>)>>,
        y_order: Option<&str>,
        decreasing: bool,
        subsamples: Option<usize>,
        subsample_size: Option<Either<usize, f64>>,
        replace: bool,
        settings: Option<HashMap<String, Py<PyAny>>>,
        seed: Option<u64>,
        n_jobs: usize,
        progress: bool,
    ) -> PyResult<Self> {
        let mut covariates_allocation = None;
        let x_parsed = parse_covariates(&X, &mut covariates_allocation)?;

        assert_safe_view_f64(&y)?;
        if y.is_empty() {
            return Err(PyValueError::new_err("y is empty, need at least some data"));
        }
        let mut responses_allocation = None;
        let y_parsed = maybe_allocate(&y, &mut responses_allocation);

        let mut observed_allocation = None;
        let maybe_observed = y_observed
            .as_ref()
            .map(|array_like| maybe_allocate(array_like, &mut observed_allocation));

        let mut weight_allocation = None;
        let maybe_weights = sample_weight
            .as_ref()
            .map(|array_like| maybe_allocate(array_like, &mut weight_allocation));

        let covariate_order = X_order
            .map(|orders| {
                CovariateGroups::parse(orders, x_parsed.dimension).map_err(|e| {
                    PyValueError::new_err(format!("covariate groups couldn't be parsed: {e}"))
                })
            })
            .transpose()?;
        // covariate_order is not used anymore if the covariate is one-dimensional

        let response_order = y_order
            .map(|name| {
                name.parse().map_err(|e| {
                    PyValueError::new_err(format!("response order name couldn't be parsed: {e}"))
                })
            })
            .transpose()?
            .unwrap_or(StochasticOrder::StochasticDominance);

        let settings =
            parse_config(settings, py).map_err(|e| PyValueError::new_err(e.to_string()))?;

        let progress_bar = progress.then(|| {
            ensure_notebook_mode_set(py);
            KdamProgress::new()
        });
        let progress: &dyn ProgressTracker = match &progress_bar {
            Some(pb) => pb,
            None => &NoProgress,
        };

        let fit = py.detach(|| {
            let config = subagging::Config::parse(
                subsamples,
                subsample_size,
                replace,
                x_parsed.n,
                seed,
                n_jobs,
            )?;
            match (x_parsed.dimension, covariate_order, settings) {
                (0, _, _) => unreachable!(),
                (1, None, Right(settings) | Both(_, settings)) => {
                    subagging::Fit::<total_order::Fit>::fit(
                        x_parsed.slice,
                        y_parsed,
                        maybe_observed,
                        maybe_weights,
                        Increasing,
                        response_order,
                        decreasing,
                        (config, settings),
                        progress,
                    )
                    .map(|fit| Fit::Total {
                        fit,
                        squeeze: x_parsed.squeeze,
                    })
                }
                (_, covariate_groups, Left(settings) | Both(settings, _)) => {
                    let order = covariate_groups
                        .unwrap_or_else(|| CovariateGroups::empty(x_parsed.dimension));
                    subagging::Fit::<partial_order::Fit>::fit(
                        x_parsed.slice,
                        y_parsed,
                        maybe_observed,
                        maybe_weights,
                        order,
                        response_order,
                        decreasing,
                        (config, settings),
                        progress,
                    )
                    .map(Fit::Partial)
                }
                _ => Err(Error::CovariateDimensionMismatch {
                    shape: if x_parsed.squeeze {
                        vec![x_parsed.n]
                    } else {
                        vec![x_parsed.n, x_parsed.dimension]
                    },
                    message: "found settings which relate to a different covariate dimension",
                }),
            }
        });

        if let Some(pb) = &progress_bar {
            pb.finish();
        }

        fit.map(|fit| IDR { inner: fit })
            .map_err(|validation_error| PyValueError::new_err(validation_error.to_string()))
    }

    /// Predict the conditional mean.
    ///
    /// Parameters
    /// ----------
    /// X : array_like, shape (...,) or (..., d)
    ///     Covariate values at which to predict. The shape must be
    ///     compatible with the training covariates: (...,) or (..., 1) for
    ///     univariate models depending on how the covariates were supplied
    ///     initially, and (..., d) for multivariate models. Supports
    ///     broadcasting over leading dimensions.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Predicted conditional means. The trailing covariate dimension is
    ///     consumed, so the output shape is the input shape without the last
    ///     axis (or unchanged for squeezed univariate input).
    fn predict<'py>(
        &self,
        py: Python<'py>,
        X: PyArrayLikeDyn<f64, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let cov_array = X.as_array();
        let array = match (&self.inner, cov_array.shape()) {
            (Fit::Total { fit, squeeze: true }, _) => Ok(cov_array.map(|&c| fit.mean(c))),
            (
                Fit::Total {
                    fit,
                    squeeze: false,
                },
                &[.., 1],
            ) => {
                let (_, out_shape) = cov_array.shape().split_last().unwrap();
                // Compute result in a flat array to avoid copying while dropping the last dim
                let result = X.as_array().iter().map(|&c| fit.mean(c)).collect();
                Ok(ArrayD::from_shape_vec(IxDyn(out_shape), result).unwrap())
            }
            (Fit::Partial(fit), shape) if shape.last() == Some(&fit.covariate_groups.dimension) => {
                let (_, out_shape) = shape.split_last().unwrap();
                let mut storage = None;
                let view = X.as_array();
                let result = view
                    .rows()
                    .into_iter()
                    .map(|row| {
                        let slice = maybe_allocate_view(&row, &mut storage);
                        fit.mean(slice)
                    })
                    .collect();
                Ok(ArrayD::from_shape_vec(IxDyn(out_shape), result).unwrap())
            }
            (Fit::Total { squeeze: false, .. }, shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., 1)",
            }),
            (Fit::Partial(_), shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., d)",
            }),
        };

        finalize(array, py)
    }

    /// Predict the full CDF at all training thresholds.
    ///
    /// Parameters
    /// ----------
    /// X : array_like, shape (...,) or (..., d)
    ///     Covariate values at which to predict. The shape must be
    ///     compatible with the training covariates: (...,) or (..., 1) for
    ///     univariate models depending on how the covariates were supplied
    ///     initially, and (..., d) for multivariate models. Supports
    ///     broadcasting over leading dimensions.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     CDF values with shape ``(*input_shape_without_cov_dim, n_thresholds)``.
    ///     The last axis corresponds to the thresholds available via the
    ///     ``thresholds`` property.
    fn cdf<'py>(
        &self,
        py: Python<'py>,
        X: PyArrayLikeDyn<f64, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let covariates_array = X.as_array();
        let in_shape = covariates_array.shape();
        let array = py.detach(|| match (&self.inner, in_shape) {
            (Fit::Total { fit, squeeze: true }, _) => {
                let mut out_shape = in_shape.to_vec();
                out_shape.push(fit.n_threshold());
                let result = covariates_array.iter().flat_map(|&c| fit.cdf(c)).collect();
                Ok(ArrayD::from_shape_vec(IxDyn(&out_shape), result).unwrap())
            }
            (
                Fit::Total {
                    fit,
                    squeeze: false,
                },
                &[.., 1],
            ) => {
                let (_, base_shape) = in_shape.split_last().unwrap();
                let mut out_shape = base_shape.to_vec();
                out_shape.push(fit.n_threshold());
                let result = covariates_array.iter().flat_map(|&c| fit.cdf(c)).collect();
                Ok(ArrayD::from_shape_vec(IxDyn(&out_shape), result).unwrap())
            }
            (Fit::Partial(fit), shape) if shape.last() == Some(&fit.covariate_groups.dimension) => {
                let (_, base_shape) = in_shape.split_last().unwrap();
                let mut out_shape = base_shape.to_vec();
                out_shape.push(fit.n_threshold());
                let mut result = Vec::with_capacity(out_shape.iter().product());
                let mut storage = None;
                for row in covariates_array.rows() {
                    let slice = maybe_allocate_view(&row, &mut storage);
                    result.extend(fit.cdf(slice))
                }
                Ok(ArrayD::from_shape_vec(IxDyn(&out_shape), result).unwrap())
            }
            (Fit::Total { squeeze: false, .. }, shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., 1)",
            }),
            (Fit::Partial(_), shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., d)",
            }),
        });

        finalize(array, py)
    }

    /// Evaluate the predicted CDF at specific response values.
    ///
    /// Unlike ``cdf``, which returns the CDF at all training thresholds,
    /// this method evaluates the CDF at arbitrary response values by
    /// interpolation.
    ///
    /// Parameters
    /// ----------
    /// X : array_like, shape (...,) or (..., d)
    ///     Covariate values at which to predict. The shape must be
    ///     compatible with the training covariates: (...,) or (..., 1) for
    ///     univariate models depending on how the covariates were supplied
    ///     initially, and (..., d) for multivariate models. Supports
    ///     broadcasting over leading dimensions. The last axis is consumed
    ///     as the covariate dimension (omitted for squeezed univariate
    ///     models).
    /// y : array_like
    ///     Response values at which to evaluate the CDF. Broadcast with the
    ///     covariate array (after the covariate dimension is consumed).
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     CDF values ``P(Y <= y | X)`` with shape determined by
    ///     broadcasting the covariate and response arrays.
    fn cdf_at<'py>(
        &self,
        py: Python<'py>,
        X: PyArrayLikeDyn<f64, AllowTypeChange>,
        y: PyArrayLikeDyn<f64, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let cov = X.as_array();
        let thr = y.as_array();

        fn execute<I: CovariateInterpolator>(
            cov_interpolated: ArrayD<I>,
            thr: ArrayD<ResponseCoordinate>,
        ) -> Result<ArrayD<f64>, Error> {
            let broadcasted = broadcast(cov_interpolated.shape(), thr.shape())?;
            let mut output = ArrayD::zeros(IxDyn(&broadcasted));
            Zip::from(&mut output)
                .and_broadcast(&cov_interpolated)
                .and_broadcast(&thr)
                .for_each(|out, c, &t| {
                    *out = c.interpolate(t);
                });
            Ok(output)
        }

        let result = py.detach(|| match (&self.inner, cov.shape()) {
            (Fit::Total { fit, squeeze: true }, _) => {
                let cov_coordinates = cov.map(|&c| fit.interpolate_covariate(c));
                let thr_coordinates = thr.map(|&t| fit.get_response_coordinate(t));
                execute(cov_coordinates, thr_coordinates)
            }
            (
                Fit::Total {
                    fit,
                    squeeze: false,
                },
                &[.., 1],
            ) => {
                let (_, without_last) = cov.shape().split_last().unwrap();
                let cov_coordinates_flat =
                    cov.iter().map(|&c| fit.interpolate_covariate(c)).collect();
                let cov_coordinates =
                    ArrayD::from_shape_vec(IxDyn(without_last), cov_coordinates_flat).unwrap();
                let thr_coordinates = thr.map(|&t| fit.get_response_coordinate(t));
                execute(cov_coordinates, thr_coordinates)
            }
            (Fit::Partial(fit), shape) if shape.last() == Some(&fit.covariate_groups.dimension) => {
                let (_, without_last) = cov.shape().split_last().unwrap();
                let mut storage = None;
                let mut workspace = PredictionWorkspace::new();
                let cov_coordinates_flat = cov
                    .rows()
                    .into_iter()
                    .map(|c| {
                        let slice = maybe_allocate_view(&c, &mut storage);
                        fit.interpolate_covariate_with_workspace(slice, &mut workspace)
                    })
                    .collect();
                let cov_coordinates =
                    ArrayD::from_shape_vec(without_last, cov_coordinates_flat).unwrap();
                let thr_coordinates = thr.map(|&t| fit.get_response_coordinate(t));
                execute(cov_coordinates, thr_coordinates)
            }
            (Fit::Total { squeeze: false, .. }, shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., 1)",
            }),
            (Fit::Partial(_), shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., d)",
            }),
        });

        finalize(result, py)
    }

    /// Evaluate the predicted CDF on a grid of sorted covariates and sorted thresholds.
    ///
    /// This is a fast path for univariate (squeezed, 1-D) models that evaluates the CDF for every
    /// combination of covariate and threshold value.
    ///
    /// Parameters
    /// ----------
    /// X : array_like, shape (m,)
    ///     Covariate values (1-D).
    /// y : array_like, shape (k,)
    ///     Response (threshold) values (1-D).
    ///
    /// Returns
    /// -------
    /// numpy.ndarray, shape (m, k)
    ///     ``out[i, j]`` is ``P(Y <= y[j] | X = X[i])``.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the model was fitted on multivariate covariates.
    fn cdf_grid<'py>(
        &self,
        py: Python<'py>,
        X: PyArrayLike1<f64, AllowTypeChange>,
        y: PyArrayLike1<f64, AllowTypeChange>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let Fit::Total { fit, squeeze: true } = &self.inner else {
            return Err(PyValueError::new_err(Error::OrderMismatch.to_string()));
        };

        let cov = X.as_array();
        let thr = y.as_array();
        let flat = fit
            .predict_grid(cov.iter().copied(), thr.iter().copied())
            .collect();

        let result = Array2::from_shape_vec([cov.len(), thr.len()], flat).unwrap();
        Ok(result.into_pyarray(py))
    }

    /// Compute quantiles of the predicted distribution.
    ///
    /// Parameters
    /// ----------
    /// X : array_like, shape (...,) or (..., d)
    ///     Covariate values at which to predict. The shape must be
    ///     compatible with the training covariates: (...,) or (..., 1) for
    ///     univariate models depending on how the covariates were supplied
    ///     initially, and (..., d) for multivariate models. Supports
    ///     broadcasting over leading dimensions. The last axis is consumed
    ///     as the covariate dimension (omitted for squeezed univariate
    ///     models).
    /// q : array_like
    ///     Quantile levels in [0, 1]. Broadcast with the covariate array
    ///     (after the covariate dimension is consumed).
    /// upper : bool, optional
    ///     If False (default), return the lower quantile (left-continuous
    ///     inverse). If True, return the upper quantile (right-continuous
    ///     inverse).
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Quantile values with shape determined by broadcasting the
    ///     covariate and quantile-level arrays.
    #[pyo3(
        signature = (
            X,
            q,
            upper=false,
        )
    )]
    fn quantile<'py>(
        &self,
        py: Python<'py>,
        X: PyArrayLikeDyn<f64, AllowTypeChange>,
        q: PyArrayLikeDyn<f64>,
        upper: bool,
    ) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let cov = X.as_array();
        let prb = q.as_array();

        fn execute(
            cov_interpolated: ArrayD<impl CovariateInterpolator>,
            prb: ArrayViewD<f64>,
            upper: bool,
            thresholds: &[f64],
        ) -> Result<ArrayD<f64>, Error> {
            let broadcasted = broadcast(cov_interpolated.shape(), prb.shape())?;
            let mut output = ArrayD::zeros(IxDyn(&broadcasted));
            Zip::from(&mut output)
                .and_broadcast(&cov_interpolated)
                .and_broadcast(&prb)
                .for_each(|out, c, &p| {
                    *out = quantile(c, p, upper, thresholds);
                });
            Ok(output)
        }

        let output = py.detach(|| match (&self.inner, cov.shape()) {
            (Fit::Total { fit, squeeze: true }, _) => {
                let cov_interpolated = cov.map(|&c| fit.interpolate_covariate(c));
                execute(cov_interpolated, prb, upper, &fit.thresholds)
            }
            (
                Fit::Total {
                    fit,
                    squeeze: false,
                },
                &[.., 1],
            ) => {
                let (_, without_last) = cov.shape().split_last().unwrap();
                let flat = cov.iter().map(|&c| fit.interpolate_covariate(c)).collect();
                let cov_interpolated = ArrayD::from_shape_vec(IxDyn(without_last), flat).unwrap();
                execute(cov_interpolated, prb, upper, &fit.thresholds)
            }
            (Fit::Partial(fit), shape) if shape.last() == Some(&fit.covariate_groups.dimension) => {
                let (_, without_last) = cov.shape().split_last().unwrap();
                let mut storage = None;
                let mut workspace = PredictionWorkspace::new();
                let flat = cov
                    .rows()
                    .into_iter()
                    .map(|c| {
                        let slice = maybe_allocate_view(&c, &mut storage);
                        fit.interpolate_covariate_with_workspace(slice, &mut workspace)
                    })
                    .collect();
                let cov_interpolated = ArrayD::from_shape_vec(IxDyn(without_last), flat).unwrap();
                execute(cov_interpolated, prb, upper, &fit.thresholds)
            }
            (Fit::Total { squeeze: false, .. }, shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., 1)",
            }),
            (Fit::Partial(_), shape) => Err(Error::CovariateDimensionMismatch {
                shape: shape.to_vec(),
                message: "expected an argument of shape (..., d)",
            }),
        });

        finalize(output, py)
    }

    /// The unique training covariates (read-only).
    ///
    /// For univariate (squeezed) models this is a 1-D array. For
    /// univariate models fitted with a 2-D input ``(n, 1)`` this is a
    /// 2-D array with one column. For multivariate models this is a 2-D
    /// array of shape ``(n_unique, d)``.
    #[getter]
    fn X<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArrayDyn<f64>> {
        let me = this.borrow();
        let array_view = match &me.inner {
            Fit::Partial(fit) => {
                let covariate_dimension = fit.covariate_groups.dimension;
                let n_covariate = fit.covariates.len() / covariate_dimension;
                ArrayView2::from_shape((n_covariate, covariate_dimension), &fit.covariates)
                    .unwrap()
                    .into_dyn()
            }
            Fit::Total { fit, squeeze: true } => ArrayView1::from(&fit.covariates).into_dyn(),
            Fit::Total {
                fit,
                squeeze: false,
            } => ArrayView2::from_shape((fit.covariates.len(), 1), &fit.covariates)
                .unwrap()
                .into_dyn(),
        };

        // SAFETY: data lives inside `this` and won't be reallocated. We set the pyclass `this` as
        // the base; Python keeps `this` alive while the NumPy array exists.
        let arr = unsafe { PyArrayDyn::borrow_from_array(&array_view, this.into_any()) };
        arr.readwrite().make_nonwriteable(); // now Python sees a read-only view
        arr
    }

    /// The unique training response thresholds (read-only).
    ///
    /// A sorted 1-D array of the distinct response values observed during
    /// training. These are the points at which ``cdf`` returns values.
    #[getter]
    fn thresholds<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f64>> {
        let me = this.borrow();
        let slice = match &me.inner {
            Fit::Partial(fit) => &fit.thresholds,
            Fit::Total { fit, .. } => &fit.thresholds,
        };
        let array_view = ArrayView1::from_shape(slice.len(), slice).unwrap();
        // SAFETY: data lives inside `this` and won't be reallocated. We set the pyclass `this` as
        // the base; Python keeps `this` alive while the NumPy array exists.
        let arr = unsafe { PyArray1::borrow_from_array(&array_view, this.into_any()) };
        arr.readwrite().make_nonwriteable(); // now Python sees a read-only view
        arr
    }

    fn __repr__(&self) -> String {
        let (order, n_covariates, covariate_dim, n_thresholds, n_subsamples, increasing) =
            self.summary();
        let direction = if increasing {
            "increasing"
        } else {
            "decreasing"
        };
        let mut parts = vec![
            format!("covariates={n_covariates}"),
            format!("dim={covariate_dim}"),
            format!("thresholds={n_thresholds}"),
            format!("order='{order}'"),
            format!("direction='{direction}'"),
        ];
        if n_subsamples > 1 {
            parts.push(format!("subsamples={n_subsamples}"));
        }
        format!("IDR({})", parts.join(", "))
    }

    fn __str__(&self) -> String {
        let (order, n_covariates, covariate_dim, n_thresholds, n_subsamples, increasing) =
            self.summary();
        let direction = if increasing {
            "increasing"
        } else {
            "decreasing"
        };
        let mut lines = vec![
            "Isotonic Distributional Regression model".to_string(),
            format!("  Order:       {order} ({direction})"),
            format!(
                "  Covariates:  {n_covariates} unique value{} ({covariate_dim}-dimensional)",
                if n_covariates == 1 { "" } else { "s" }
            ),
            format!(
                "  Thresholds:  {n_thresholds} unique value{}",
                if n_thresholds == 1 { "" } else { "s" }
            ),
        ];
        if n_subsamples > 1 {
            lines.push(format!("  Subsamples:  {n_subsamples}"));
        }
        lines.join("\n")
    }

    // Pickle: return (callable, args)
    fn __reduce_ex__<'py>(&self, py: Python<'py>, _protocol: u8) -> PyResult<Bound<'py, PyTuple>> {
        let buf = bincode::serde::encode_to_vec(&self.inner, bincode::config::standard())
            .map_err(|e| PyException::new_err(e.to_string()))?;
        let bytes = PyBytes::new(py, &buf);
        let ctor = py.get_type::<Self>().getattr("_from_bytes")?;
        let args = PyTuple::new(py, [bytes])?;
        // Return the 2‑tuple (callable, args) as a Python tuple
        PyTuple::new(py, [ctor, args.into_any()])
    }

    #[classmethod]
    fn _from_bytes(_cls: &Bound<'_, PyType>, b: &[u8]) -> PyResult<Self> {
        let (inner, _) = bincode::serde::decode_from_slice(b, bincode::config::standard())
            .map_err(|e| PyException::new_err(e.to_string()))?;
        Ok(IDR { inner })
    }

    /// Create an IDR model from pre-computed CDFs, covariates, and thresholds.
    ///
    /// Parameters
    /// ----------
    /// cdfs : array_like, shape (n_covariates, n_thresholds)
    ///     Matrix of CDF values (one row per unique covariate).
    /// X : array_like, shape (n_covariates,) or (n_covariates, d)
    ///     Unique covariate values in strictly increasing lexicographic order.
    /// y : array_like, shape (n_thresholds,)
    ///     Unique thresholds in strictly increasing order.
    /// global_cdf : array_like, shape (n_thresholds,), optional
    ///     For multivariate covariates, the CDF prediction for covariates
    ///     that are incomparable to every training covariate. When omitted,
    ///     approximated by the unweighted average of the per-covariate CDFs.
    ///
    /// Returns
    /// -------
    /// IDR
    ///     A fitted IDR model.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If dimensions are inconsistent, thresholds are not strictly
    ///     increasing, or covariates are not in strictly increasing
    ///     lexicographic order.
    #[classmethod]
    #[pyo3(signature = (cdfs, X, y, global_cdf=None))]
    fn from_cdfs(
        _cls: &Bound<'_, PyType>,
        cdfs: PyArrayLike2<f64, AllowTypeChange>,
        X: PyArrayLikeDyn<f64, AllowTypeChange>,
        y: PyArrayLike1<f64, AllowTypeChange>,
        global_cdf: Option<PyArrayLike1<f64, AllowTypeChange>>,
    ) -> PyResult<Self> {
        let mut covariates_allocation = None;
        let covariates = parse_covariates(&X, &mut covariates_allocation)?;

        let cdfs_array = cdfs.as_array();
        let thresholds_array = y.as_array();

        let n_cdfs_rows = cdfs_array.shape()[0];
        let n_cdfs_cols = cdfs_array.shape()[1];
        let n_thresh = thresholds_array.len();

        // Dimension checks
        if n_cdfs_rows != covariates.n {
            let n = covariates.n;
            return Err(PyValueError::new_err(format!(
                "Number of CDF rows ({n_cdfs_rows}) must match number of covariates ({n})"
            )));
        }
        if n_cdfs_cols != n_thresh {
            return Err(PyValueError::new_err(format!(
                "Number of CDF columns ({n_cdfs_cols}) must match number of thresholds ({n_thresh})"
            )));
        }

        let thresholds_vec: Vec<f64> = thresholds_array.iter().copied().collect();
        if !thresholds_vec.windows(2).all(|w| w[0] < w[1]) {
            return Err(PyValueError::new_err(
                "Thresholds must be sorted in strictly increasing order",
            ));
        }

        // Check covariates are sorted in strictly increasing lexicographic order
        for (left, right) in covariates
            .slice
            .chunks_exact(covariates.dimension)
            .zip(covariates.slice.chunks_exact(covariates.dimension).skip(1))
        {
            if lexicographic_cmp(left, right) != Ordering::Less {
                return Err(PyValueError::new_err(
                    "Covariates must be sorted in strictly increasing lexicographic order and contain no duplicates",
                ));
            }
        }

        // Flatten CDFs in row-major order (= covariate-major, which is the internal layout)
        let cdfs_flat: Vec<f64> = cdfs_array.iter().copied().collect();

        let covariates_owned = covariates.slice.to_vec();

        match covariates.dimension {
            0 => unreachable!(),
            1 => {
                let subagging_covariates = covariates_owned.clone();
                let inner_fit = total_order::Fit {
                    increasing: true,
                    cdfs: cdfs_flat,
                    covariates: covariates_owned,
                    thresholds: thresholds_vec.clone(),
                    quality_indicators: total_order::QualityIndicators { epsilon: f64::NAN },
                };
                let subagging_fit = subagging::Fit::from_parts(
                    inner_fit,
                    subagging_covariates,
                    thresholds_vec,
                    Increasing,
                );
                Ok(IDR {
                    inner: Fit::Total {
                        fit: subagging_fit,
                        squeeze: covariates.squeeze,
                    },
                })
            }
            _ => {
                let covariate_groups = CovariateGroups::empty(covariates.dimension);

                let (ordering_info, edges_covariates) = if covariates.n > 0 {
                    let edges = derive_transitive_reduction(
                        &covariates_owned,
                        covariates.n,
                        covariates.dimension,
                    );
                    (
                        OrderingInfo::from_edges(edges, covariates.n),
                        covariates_owned.clone(),
                    )
                } else {
                    (OrderingInfo::empty(), Vec::new())
                };

                let global_cdf_vec = if let Some(ref gcdf) = global_cdf {
                    let gcdf_array = gcdf.as_array();
                    if gcdf_array.len() != n_thresh {
                        return Err(PyValueError::new_err(format!(
                            "global_cdf length ({}) must match number of thresholds ({n_thresh})",
                            gcdf_array.len(),
                        )));
                    }
                    gcdf_array.iter().copied().collect()
                } else if covariates.n > 0 {
                    // Unweighted average of per-covariate CDFs; a rough approximation
                    // since the true global CDF should be weighted by the number of
                    // observations at each covariate.
                    let mut avg = vec![0.0; n_thresh];
                    for i in 0..covariates.n {
                        for j in 0..n_thresh {
                            avg[j] += cdfs_flat[i * n_thresh + j];
                        }
                    }
                    for v in &mut avg {
                        *v /= covariates.n as f64;
                    }
                    avg
                } else {
                    Vec::new()
                };

                let inner_fit = partial_order::Fit {
                    increasing: true,
                    cdfs: cdfs_flat,
                    global_cdf: global_cdf_vec,
                    covariate_groups: covariate_groups.clone(),
                    covariates: covariates_owned,
                    ordering_info,
                    thresholds: thresholds_vec.clone(),
                    quality_indicators: QualityIndicators {
                        precision: f64::NAN,
                        convergence_fraction: f64::NAN,
                    },
                };
                let subagging_fit = subagging::Fit::from_parts(
                    inner_fit,
                    edges_covariates,
                    thresholds_vec,
                    covariate_groups,
                );
                Ok(IDR {
                    inner: Fit::Partial(subagging_fit),
                })
            }
        }
    }
}

/// kdam-backed progress tracker.
///
/// `Mutex<Bar>` provides `Send + Sync`. kdam internally rate-limits redraws so calls to
/// `update(1)` are dominated by an integer increment; lock contention is negligible compared
/// to the per-threshold work that wraps each call.
struct KdamProgress {
    bar: Mutex<Bar>,
}

impl KdamProgress {
    fn new() -> Self {
        Self {
            bar: Mutex::new(tqdm!(total = 0, force_refresh = true, desc = "Fitting IDR")),
        }
    }
}

impl ProgressTracker for KdamProgress {
    fn set_total(&self, n: usize) {
        if let Ok(mut bar) = self.bar.lock() {
            bar.total = n;
        }
    }
    /// Should not be called too frequently.
    fn increment(&self) {
        if let Ok(mut bar) = self.bar.lock() {
            let _ = bar.update(1);
        }
    }

    fn finish(&self) {
        if let Ok(mut bar) = self.bar.lock() {
            let _ = bar.refresh();
        }
    }
}

/// Detect a Jupyter / IPython kernel and switch kdam into widget rendering mode the first
/// time we observe one in the current process. Without this, kdam falls back to `\r`-based
/// stderr output which renders as a wall of duplicate lines in notebooks.
fn ensure_notebook_mode_set(py: Python) {
    static DONE: OnceLock<()> = OnceLock::new();
    DONE.get_or_init(|| {
        let in_jupyter = py
            .import("sys")
            .and_then(|sys| sys.getattr("modules"))
            .and_then(|m| m.contains("ipykernel"))
            .unwrap_or(false);
        if in_jupyter {
            kdam::set_notebook(true);
        }
    });
}

fn finalize(
    result: Result<Array<f64, IxDyn>, Error>,
    py: Python,
) -> PyResult<Bound<PyArrayDyn<f64>>> {
    match result {
        Ok(array) => Ok(array.into_pyarray(py)),
        Err(e) => Err(PyValueError::new_err(e.to_string())),
    }
}

fn broadcast(covariate: &[usize], other: &[usize]) -> Result<Vec<usize>, Error> {
    broadcast_shapes([covariate, other]).ok_or_else(|| Error::ShapeMismatch {
        covariate_shape: covariate.to_vec(),
        other_shape: other.to_vec(),
    })
}

fn parse_covariates<'a>(
    covariates: &'a PyArrayLikeDyn<f64, AllowTypeChange>,
    storage: &'a mut Option<Vec<f64>>,
) -> PyResult<Covariates<'a>> {
    assert_safe_view_f64(covariates)?;

    match covariates.shape() {
        &[0] => Err(PyValueError::new_err(
            "covariates: Expected a non-empty array",
        )),
        &[n] => Ok(Covariates {
            slice: maybe_allocate(covariates, storage),
            n,
            dimension: 1,
            squeeze: true,
        }),
        &[_, 0] => Err(PyValueError::new_err(
            "covariates: Last dimension has size 0, expected a non-empty array instead",
        )),
        &[0, _] => Err(PyValueError::new_err(
            "covariates: Expected a non-empty array",
        )),
        &[n, dimension] => Ok(Covariates {
            slice: maybe_allocate(covariates, storage),
            n,
            dimension,
            squeeze: false,
        }),
        shape => Err(PyValueError::new_err(format!(
            "covariates: Expected a 1D or 2D array, got {}D array instead",
            shape.len(),
        ))),
    }
}

/// Represents a 2d array.
///
/// TODO: Replace with something from an existing library?
struct Covariates<'a> {
    slice: &'a [f64],
    n: usize,
    dimension: usize,
    squeeze: bool,
}

fn parse_config(
    maybe_config: Option<HashMap<String, Py<PyAny>>>,
    py: Python,
) -> Result<EitherOrBoth<partial_order::Config, total_order::Config>, Error> {
    let mut partial_order_config = partial_order::Config::default();
    let mut total_order_config = total_order::Config::default();

    let Some(mut config) = maybe_config else {
        return Ok(Both(partial_order_config, total_order_config));
    };

    // Track whether any of the configs is modified
    let mut any_partial_order_options = false;
    let mut any_total_order_options = false;

    // No fields are shared between them

    // Partial order config fields
    if let Some(map) = config.remove("osqp_settings") {
        any_partial_order_options = true;
        let mut settings = mem::take(&mut partial_order_config.osqp_settings);
        let settings_map = map.extract::<HashMap<String, Py<PyAny>>>(py).map_err(|_| {
            Error::ConfigParseError("osqp_settings should be a dictionary with string keys")
        })?;

        for (key, value) in &settings_map {
            match key.as_str() {
                "verbose" => {
                    let v = value
                        .extract::<bool>(py)
                        .map_err(|_| Error::ConfigParseError("verbose should be a boolean"))?;
                    settings = settings.verbose(v);
                }
                "eps_abs" => {
                    let v = value
                        .extract::<f64>(py)
                        .map_err(|_| Error::ConfigParseError("eps_abs should be a float"))?;
                    settings = settings.eps_abs(v);
                }
                "eps_rel" => {
                    let v = value
                        .extract::<f64>(py)
                        .map_err(|_| Error::ConfigParseError("eps_rel should be a float"))?;
                    settings = settings.eps_rel(v);
                }
                "max_iter" => {
                    let v = value.extract::<u32>(py).map_err(|_| {
                        Error::ConfigParseError("max_iter should be a positive integer")
                    })?;
                    settings = settings.max_iter(v);
                }
                _ => {
                    return Err(Error::ConfigParseError(
                        "unknown osqp_settings key, they need to be added individually in the isodistrreg rust package",
                    ));
                }
            }
        }

        partial_order_config.osqp_settings = settings;
    }

    // Total order config fields
    if let Some(epsilon) = config.remove("epsilon") {
        any_total_order_options = true;
        let as_f64 = epsilon
            .extract::<f64>(py)
            .map_err(|_| {
                Error::ConfigParseError("could not convert provided epsilon to a positive float")
            })
            .and_then(|e| {
                if e > 0.0 {
                    Ok(e)
                } else {
                    Err(Error::ConfigParseError(
                        "epsilon should be strictly positive",
                    ))
                }
            })?;
        total_order_config.epsilon = as_f64;
    }

    match (any_partial_order_options, any_total_order_options) {
        (false, false) => Ok(Both(partial_order_config, total_order_config)),
        (true, false) => Ok(Left(partial_order_config)),
        (false, true) => Ok(Right(total_order_config)),
        (true, true) => Err(Error::ConfigParseError(
            "mixing partial and total order settings",
        )),
    }
}

/// Compute an isotonic regression for the mean.
///
/// Fits values that minimize a weighted squared error subject to a monotonicity
/// constraint over the **last** broadcast axis. Outer broadcast axes describe
/// independent regressions, and the broadcast shape determines the output shape.
///
/// One of four schemes is chosen automatically:
///
/// 1. **Indexed total order** — `covariates` and `constraints` both ``None``.
///    Each cell along the last broadcast axis is ordered by its index.
/// 2. **1-D covariate** — `covariates` broadcasts together with `responses`
///    and `weights`. Each cell is one observation. Equal covariate values are
///    pooled, so the output's last axis can be shorter than the input's.
/// 3. **Multidimensional covariate** — `covariates` has a trailing axis of
///    length ``d > 1`` that does not broadcast with `responses` (or `weights`).
///    Each ``d``-row is one observation, compared componentwise.
/// 4. **Custom partial order** — `constraints` is an ``(m, 2)`` array of index
///    pairs ``(i, j)`` meaning "observation ``i`` precedes ``j``". The same
///    edge set is applied to every regression along the last broadcast axis.
///
/// `covariates` and `constraints` are mutually exclusive. In cases 2 and 3,
/// outputs are returned in the order of the (sorted, lexicographic for case 3)
/// unique covariate values.
///
/// Parameters
/// ----------
/// y : array_like
///     Response values; the last broadcast axis is the regression axis.
/// X : array_like, optional
///     Order-defining covariates. A trailing axis of length ``d > 1`` triggers
///     the multidimensional case.
/// sample_weight : array_like, optional
///     Non-negative observation weights. Defaults to uniform weights.
/// decreasing : bool, default False
///     Fit a monotone-decreasing regression instead of increasing.
/// constraints : array_like of shape (m, 2), optional
///     Partial-order edges as ``(i, j)`` index pairs.
///
/// Returns
/// -------
/// numpy.ndarray of float64
///     Fitted values, broadcast across outer axes.
///
/// Raises
/// ------
/// ValueError
///     If shapes cannot be broadcast, if both `covariates` and `constraints`
///     are given, or if independent regressions over the outer broadcast axes
///     produce outputs of different lengths (see last example).
///
/// Examples
/// --------
/// **Case 1** — indexed total order. Pools the non-monotone sequence to its
/// running mean:
///
///     >>> isotonic_regression([4.0, 3.0, 2.0])
///     array([3., 3., 3.])
///
/// **Case 2** — 1-D covariate. Observations are paired with covariate values;
/// equal covariates pool, and the output is in sorted-covariate order:
///
///     >>> isotonic_regression([4.0, 3.0, 2.0], covariates=[2.0, 1.0, 1.0])
///     array([2.5, 4. ])
///
/// **Case 3** — multidimensional covariate. Each row of `covariates` is a
/// ``d``-dim point, compared componentwise; outputs follow lexicographic order
/// of unique rows:
///
///     >>> covariates = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
///     >>> isotonic_regression([1.0, 2.0, 3.0, 4.0], covariates=covariates)
///     array([1., 3., 2., 4.])
///
/// **Case 4** — explicit partial order. Here, the single edge ``1 -> 0`` forces
/// the second observation to precede the first, pooling them:
///
///     >>> isotonic_regression([1.0, 2.0], constraints=[(1, 0)])
///     array([1.5, 1.5])
///
/// **Broadcasting** — run a stack of regressions in parallel. With responses
/// of shape ``(2, 3)``, two independent length-3 regressions run:
///
///     >>> import numpy as np
///     >>> isotonic_regression(np.array([[4.0, 3.0, 2.0], [1.0, 2.0, 3.0]]))
///     array([[3., 3., 3.],
///            [1., 2., 3.]])
///
/// Outer axes can come from any input; the result takes the broadcast shape:
///
///     >>> isotonic_regression(
///     ...     responses=np.arange(12).reshape(3, 4),
///     ...     weights=np.arange(1, 13).reshape(1, 1, 3, 4),
///     ... ).shape
///     (1, 1, 3, 4)
///
/// **Output-length mismatch** — covariate pooling is per-regression, so a
/// regression whose covariates contain duplicates produces a shorter output
/// than one with all unique values. When the resulting outputs would not
/// stack into a rectangular array, the call is rejected:
///
///     >>> responses = np.arange(12).reshape(3, 4)
///     >>> covariates = np.stack([np.arange(4), np.ones(4), np.arange(4)])
///     >>> isotonic_regression(responses, covariates=covariates)
///     Traceback (most recent call last):
///         ...
///     ValueError: not all regression results have the same length
#[allow(non_snake_case)]
#[pyfunction(
    signature = (
        y,
        X=None,
        sample_weight=None,
        decreasing=false,
        constraints=None,
    )
)]
fn isotonic_regression<'py>(
    py: Python<'py>,
    y: PyArrayLikeDyn<'py, f64, AllowTypeChange>,
    X: Option<PyArrayLikeDyn<'py, f64, AllowTypeChange>>,
    sample_weight: Option<PyArrayLikeDyn<'py, f64, AllowTypeChange>>,
    decreasing: bool,
    constraints: Option<PyArrayLike2<'py, usize, AllowTypeChange>>,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    if X.is_some() && constraints.is_some() {
        return Err(PyValueError::new_err(
            "supply either totally ordered covariates or partial order constraints, not both",
        ));
    }

    let responses = y.as_array();
    let unit_weight = ArrayD::ones(IxDyn(&[1]));
    let weights = sample_weight
        .as_ref()
        .map(|w| w.as_array())
        .unwrap_or_else(|| unit_weight.view());

    if let Some(covariates_like) = X.as_ref() {
        let covariates = covariates_like.as_array();
        // First try to broadcast all three arrays together (covariate treated as 1-D).
        // On failure, treat the covariate's last axis as the covariate dimension.
        if let Some(shape) =
            broadcast_shapes([covariates.shape(), responses.shape(), weights.shape()])
        {
            isotonic_1d_covariate(py, covariates, responses, weights, shape, decreasing)
        } else {
            isotonic_multidim_covariate(py, covariates, responses, weights, decreasing)
        }
    } else if let Some(constraints) = constraints.as_ref() {
        isotonic_constrained(py, responses, weights, constraints.as_array(), decreasing)
    } else {
        isotonic_indexed(py, responses, weights, decreasing)
    }
}

fn broadcast_err(names: &str, shapes: &[&[usize]]) -> PyErr {
    PyValueError::new_err(format!(
        "{names} could not be broadcast together with shapes {shapes:?}"
    ))
}

/// No covariate, no constraint: each cell along the last broadcast axis is indexed in order.
/// Output shape equals the broadcast shape, so we preallocate and fill in place.
fn isotonic_indexed<'py>(
    py: Python<'py>,
    y: ArrayViewD<'_, f64>,
    sampel_weight: ArrayViewD<'_, f64>,
    decreasing: bool,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    fn run<D: Direction>(
        output: &mut ArrayD<f64>,
        responses: ArrayViewD<'_, f64>,
        weights: ArrayViewD<'_, f64>,
        axis: Axis,
    ) {
        Zip::from(output.lanes_mut(axis))
            .and(responses.lanes(axis))
            .and(weights.lanes(axis))
            .for_each(|o, r, w| {
                write_lane(
                    o,
                    total_order::tonic_regression_pre_sorted::<D>(
                        r.iter().copied().zip(w.iter().copied()),
                    ),
                );
            });
    }

    let shape = broadcast_shapes([y.shape(), sampel_weight.shape()]).ok_or_else(|| {
        broadcast_err("responses and weights", &[y.shape(), sampel_weight.shape()])
    })?;
    let axis = Axis(shape.len() - 1);
    let mut output = ArrayD::zeros(IxDyn(&shape));
    let responses = y.broadcast(IxDyn(&shape)).unwrap();
    let weights = sampel_weight.broadcast(IxDyn(&shape)).unwrap();

    match decreasing {
        false => run::<Increasing>(&mut output, responses, weights, axis),
        true => run::<Decreasing>(&mut output, responses, weights, axis),
    }
    Ok(output.into_pyarray(py))
}

/// User-supplied partial-order edges on the last broadcast axis. Responses and weights broadcast;
/// the same edge set is applied to each lane. Output shape equals the broadcast shape.
fn isotonic_constrained<'py>(
    py: Python<'py>,
    y: ArrayViewD<'_, f64>,
    weight: ArrayViewD<'_, f64>,
    constraint: ArrayView2<'_, usize>,
    decreasing: bool,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    fn run<D: Direction>(
        output: &mut ArrayD<f64>,
        responses: ArrayViewD<'_, f64>,
        weights: ArrayViewD<'_, f64>,
        edges: &[(usize, usize)],
        axis: Axis,
    ) {
        Zip::from(output.lanes_mut(axis))
            .and(responses.lanes(axis))
            .and(weights.lanes(axis))
            .for_each(|o, r, w| {
                let obs = r.iter().zip(w.iter()).map(|(&r, &w)| Observation {
                    x: (),
                    y: r,
                    observed: (),
                    weight: w,
                });
                write_lane(
                    o,
                    partial_order::tonic_regression_pre_sorted::<D, _, _>(
                        obs,
                        edges,
                        &Average::new(),
                    ),
                );
            });
    }

    if constraint.shape()[1] != 2 {
        return Err(PyValueError::new_err(
            "constraints must be of dimension (m, 2)",
        ));
    }
    let shape = broadcast_shapes([y.shape(), weight.shape()])
        .ok_or_else(|| broadcast_err("responses and weights", &[y.shape(), weight.shape()]))?;
    let axis = Axis(shape.len().max(1) - 1);
    let mut output = ArrayD::zeros(IxDyn(&shape));
    let responses = y.broadcast(IxDyn(&shape)).unwrap();
    let weights = weight.broadcast(IxDyn(&shape)).unwrap();
    let edges: Vec<_> = constraint
        .outer_iter()
        .map(|row| (row[0], row[1]))
        .collect();

    match decreasing {
        false => run::<Increasing>(&mut output, responses, weights, &edges, axis),
        true => run::<Decreasing>(&mut output, responses, weights, &edges, axis),
    }
    Ok(output.into_pyarray(py))
}

/// Covariate broadcasts alongside responses and weights: each cell along the last broadcast
/// axis is one `(covariate, response, weight)` observation, and that axis is regressed over.
/// Duplicate covariate values collapse, so the output axis can be shorter than the input.
fn isotonic_1d_covariate<'py>(
    py: Python<'py>,
    x: ArrayViewD<'_, f64>,
    y: ArrayViewD<'_, f64>,
    weight: ArrayViewD<'_, f64>,
    shape: Vec<usize>,
    decreasing: bool,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    fn run<D: Direction>(
        x: ArrayViewD<'_, f64>,
        y: ArrayViewD<'_, f64>,
        weight: ArrayViewD<'_, f64>,
        axis: Axis,
    ) -> PyResult<(Vec<f64>, usize)> {
        let mut flat = Vec::new();
        let mut output_length: Option<usize> = None;
        let mut err: Option<PyErr> = None;
        Zip::from(x.lanes(axis))
            .and(y.lanes(axis))
            .and(weight.lanes(axis))
            .for_each(|c, r, w| {
                if err.is_some() {
                    return;
                }
                let iter = total_order::tonic_regression::<D>(izip!(
                    c.iter().copied(),
                    r.iter().copied(),
                    w.iter().copied(),
                ));
                if let Err(e) = reconcile_length(&mut output_length, iter.len()) {
                    err = Some(e);
                    return;
                }
                flat.extend(iter);
            });
        if let Some(e) = err {
            return Err(e);
        }
        Ok((flat, output_length.unwrap_or(0)))
    }

    let axis = Axis(shape.len() - 1);
    let covariates = x.broadcast(IxDyn(&shape)).unwrap();
    let responses = y.broadcast(IxDyn(&shape)).unwrap();
    let weights = weight.broadcast(IxDyn(&shape)).unwrap();

    let (flat, output_length) = match decreasing {
        false => run::<Increasing>(covariates, responses, weights, axis)?,
        true => run::<Decreasing>(covariates, responses, weights, axis)?,
    };

    let mut result_shape = shape;
    *result_shape.last_mut().unwrap() = output_length;
    Ok(ArrayD::from_shape_vec(IxDyn(&result_shape), flat)
        .unwrap()
        .into_pyarray(py))
}

/// Covariate has a trailing axis of length `d > 1`: each `(..., d)` row is one multidimensional
/// covariate, compared componentwise. Responses and weights broadcast with everything except
/// that last axis. The output axis is the number of unique covariate rows per regression.
fn isotonic_multidim_covariate<'py>(
    py: Python<'py>,
    x: ArrayViewD<'_, f64>,
    y: ArrayViewD<'_, f64>,
    weight: ArrayViewD<'_, f64>,
    decreasing: bool,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    let c_shape = x.shape();
    let covariate_dimension = *c_shape.last().unwrap_or(&1);
    let broadcast_err = || {
        broadcast_err(
            "covariates, responses, and weights",
            &[c_shape, y.shape(), weight.shape()],
        )
    };
    if c_shape.len() < 2 || covariate_dimension <= 1 {
        return Err(broadcast_err());
    }
    let without_last = c_shape.split_last().unwrap().1;
    let outer_shape =
        broadcast_shapes([without_last, y.shape(), weight.shape()]).ok_or_else(broadcast_err)?;

    // Broadcast everything to at least 1-D outer space so `exact_chunks` has room for the
    // regression axes.
    let mut scalar_shape = outer_shape.clone();
    if scalar_shape.is_empty() {
        scalar_shape.push(1);
    }
    let mut cov_shape = scalar_shape.clone();
    cov_shape.push(covariate_dimension);

    // Align the rank of responses and weights with covariates (which has an extra `d` axis) so
    // that `Zip` can iterate them in lockstep.
    let covariates = x.broadcast(IxDyn(&cov_shape)).unwrap();
    let responses = y
        .broadcast(IxDyn(&scalar_shape))
        .unwrap()
        .insert_axis(Axis(scalar_shape.len()));
    let weights = weight
        .broadcast(IxDyn(&scalar_shape))
        .unwrap()
        .insert_axis(Axis(scalar_shape.len()));

    // Each chunk is one regression; outer axes get size-1 placeholders.
    let mut cov_chunk = cov_shape;
    let mut scalar_chunk = scalar_shape;
    scalar_chunk.push(1);
    let cov_outer = cov_chunk.len() - 2;
    cov_chunk[..cov_outer].fill(1);
    scalar_chunk[..cov_outer].fill(1);

    let mut flat = Vec::new();
    let mut output_length: Option<usize> = None;
    let mut err: Option<PyErr> = None;
    Zip::from(covariates.exact_chunks(IxDyn(&cov_chunk)))
        .and(responses.exact_chunks(IxDyn(&scalar_chunk)))
        .and(weights.exact_chunks(IxDyn(&scalar_chunk)))
        .for_each(|c, r, w| {
            if err.is_some() {
                return;
            }
            let (mut c_owned, mut r_owned, mut w_owned) = (None, None, None);
            let result = partial_order_regression(
                get_slice(&c, &mut c_owned),
                get_slice(&r, &mut r_owned),
                get_slice(&w, &mut w_owned),
                covariate_dimension,
                decreasing,
            );
            if let Err(e) = reconcile_length(&mut output_length, result.len()) {
                err = Some(e);
                return;
            }
            flat.extend(result);
        });
    if let Some(e) = err {
        return Err(e);
    }

    let mut result_shape = outer_shape;
    result_shape.pop();
    result_shape.push(output_length.unwrap_or(0));
    Ok(ArrayD::from_shape_vec(IxDyn(&result_shape), flat)
        .unwrap()
        .into_pyarray(py))
}

/// Deduplicate a multidimensional covariate and run a partial-order regression on the result.
fn partial_order_regression(
    x: &[f64],
    y: &[f64],
    weight: &[f64],
    dimension: usize,
    decreasing: bool,
) -> Vec<f64> {
    let context = preprocess_uncensored(x, y, weight, &CovariateGroups::empty(dimension));
    let n_unique = context.x.len() / dimension;
    let edges = derive_transitive_reduction(&context.x, n_unique, dimension);

    let mut response_sums = vec![0.0; n_unique];
    for (r, w, i) in izip!(&context.y, &context.weight, &context.x_indices,) {
        response_sums[*i] += w * r;
    }
    debug_assert!(context.x_weight.iter().all(|&w| w > 0.0));
    let obs = izip!(response_sums, &context.x_weight).map(|(sum, &weight)| Observation {
        x: (),
        y: sum / weight,
        observed: (),
        weight,
    });

    match decreasing {
        false => partial_order::tonic_regression_pre_sorted::<Increasing, _, _>(
            obs,
            &edges,
            &Average::new(),
        ),
        true => partial_order::tonic_regression_pre_sorted::<Decreasing, _, _>(
            obs,
            &edges,
            &Average::new(),
        ),
    }
}

/// Write the values yielded by `iter` into `out` in order; the iterator must yield at least
/// `out.len()` elements.
fn write_lane<I: IntoIterator<Item = f64>>(mut out: ArrayViewMut1<'_, f64>, iter: I) {
    let expected = out.len();
    let mut written = 0;
    for (slot, value) in out.iter_mut().zip(iter) {
        *slot = value;
        written += 1;
    }
    assert_eq!(
        written, expected,
        "write_lane: iterator yielded {written} elements, lane has {expected}"
    );
}

fn reconcile_length(output_length: &mut Option<usize>, this: usize) -> PyResult<()> {
    match *output_length {
        None => {
            *output_length = Some(this);
            Ok(())
        }
        Some(previous) if previous == this => Ok(()),
        Some(_) => Err(PyValueError::new_err(
            "not all regression results have the same length",
        )),
    }
}

fn assert_safe_view_f64<I>(array: &Bound<PyArray<f64, I>>) -> PyResult<()> {
    // 1) Endianness
    match array.as_untyped().dtype().is_native_byteorder() {
        Some(true) | None => {} // '=' or not applicable
        Some(false) => {
            return Err(PyValueError::new_err(
                "float64 array is not native-endian (byte-swapped)",
            ));
        }
    }

    let align = align_of::<f64>(); // 8 on all modern platforms

    // 2) Data pointer alignment
    let data_ptr = array.data() as usize; // PyArrayMethods::data()
    if !data_ptr.is_multiple_of(align) {
        return Err(PyValueError::new_err(format!(
            "unaligned data pointer for f64: ptr={:#x}, align={}",
            data_ptr, align
        )));
    }

    // 3) Strides alignment (bytes, may be negative)
    for (axis, &s) in array.as_untyped().strides().iter().enumerate() {
        if s.unsigned_abs() % align != 0 {
            return Err(PyValueError::new_err(format!(
                "unaligned stride on axis {}: stride={} bytes, required multiple of {}",
                axis, s, align
            )));
        }
    }

    Ok(())
}

/// Compute the NumPy-style broadcast shape of several shapes.
///
/// Returns `None` if the shapes are not broadcast-compatible.
fn broadcast_shapes<const M: usize>(shapes: [&[usize]; M]) -> Option<Vec<usize>> {
    if shapes.is_empty() {
        return Some(Vec::new());
    }

    fn broadcast_two(a: &[usize], b: &[usize]) -> Option<Vec<usize>> {
        let ndim = a.len().max(b.len());
        let mut shape = Vec::with_capacity(ndim);
        for i in 0..ndim {
            let da = if i + a.len() >= ndim {
                a[i + a.len() - ndim]
            } else {
                1
            };
            let db = if i + b.len() >= ndim {
                b[i + b.len() - ndim]
            } else {
                1
            };
            match (da, db) {
                (x, y) if x == y => shape.push(x),
                (1, y) => shape.push(y),
                (x, 1) => shape.push(x),
                _ => return None,
            }
        }
        Some(shape)
    }

    let mut result = shapes[0].to_vec();
    for &shape in &shapes[1..] {
        result = broadcast_two(&result, shape)?;
    }
    Some(result)
}

/// Time-axis element types accepted by [`kaplan_meier`].
///
/// Abstracts the total-order comparison used by the algorithm: floats use
/// [`f64::total_cmp`] / [`f32::total_cmp`] (NaN-safe), integers use [`Ord::cmp`].
trait TimeValue: Element + Copy + PartialEq + 'static {
    fn total_cmp(&self, other: &Self) -> Ordering;
}

macro_rules! impl_time_value_float {
    ($($t:ty),+ $(,)?) => {$(
        impl TimeValue for $t {
            #[inline]
            fn total_cmp(&self, other: &Self) -> Ordering {
                <$t>::total_cmp(self, other)
            }
        }
    )+};
}

macro_rules! impl_time_value_int {
    ($($t:ty),+ $(,)?) => {$(
        impl TimeValue for $t {
            #[inline]
            fn total_cmp(&self, other: &Self) -> Ordering {
                Ord::cmp(self, other)
            }
        }
    )+};
}

impl_time_value_float!(f32, f64);
impl_time_value_int!(i8, i16, i32, i64, u8, u16, u32, u64);

fn kaplan_meier_jumps<T: TimeValue>(
    y: ArrayView1<'_, T>,
    y_observed: ArrayView1<'_, bool>,
    weight: Option<ArrayView1<'_, f64>>,
) -> Vec<(T, f64)> {
    let n = y.len();
    match weight {
        Some(weights) => {
            let order = argsort_unstable_by::<Increasing, _>(
                // Sort event times before censoring times; sort zero-weight
                // observations to the end so they don't perturb the estimator.
                |a, b| match (weights[a] == 0.0, weights[b] == 0.0) {
                    (true, true) => Ordering::Equal,
                    (true, false) => Ordering::Greater,
                    (false, true) => Ordering::Less,
                    (false, false) => y[a]
                        .total_cmp(&y[b])
                        .then(y_observed[a].cmp(&y_observed[b]).reverse()),
                },
                n,
            );
            order
                .chunk_by(|&i, &j| y[i] == y[j] && y_observed[i] == y_observed[j])
                .scan((1.0, weights.sum()), |(s, total_weight), group| {
                    let group_weight: f64 = group.iter().map(|&i| weights[i]).sum();
                    let head = group[0];
                    let time = y[head];
                    let event = y_observed[head];

                    if *total_weight <= 0.0 {
                        return None;
                    }

                    let jump = event.then(|| {
                        *s *= (1.0 - group_weight / *total_weight).clamp(0.0, 1.0);
                        (time, *s)
                    });
                    *total_weight -= group_weight;
                    Some(jump)
                })
                .flatten()
                .collect()
        }
        None => {
            let order = argsort_unstable_by::<Increasing, _>(
                |a, b| {
                    y[a].total_cmp(&y[b])
                        .then(y_observed[a].cmp(&y_observed[b]).reverse())
                },
                n,
            );
            order
                .chunk_by(|&i, &j| y[i] == y[j] && y_observed[i] == y_observed[j])
                .scan((1.0, n), |(s, total_count), group| {
                    let group_count = group.len();
                    let head = group[0];
                    let time = y[head];
                    let event = y_observed[head];

                    let jump = event.then(|| {
                        *s *= (*total_count - group_count) as f64 / *total_count as f64;
                        (time, *s)
                    });
                    *total_count -= group_count;
                    Some(jump)
                })
                .flatten()
                .collect()
        }
    }
}

#[allow(clippy::type_complexity)]
#[pyfunction(
    signature = (
        y,
        y_observed,
        weight=None,
    )
)]
fn kaplan_meier<'py>(
    py: Python<'py>,
    y: &Bound<'py, PyAny>,
    y_observed: PyArrayLike1<'py, bool, AllowTypeChange>,
    weight: Option<PyArrayLike1<'py, f64, AllowTypeChange>>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyArray1<f64>>)> {
    // Normalize the input into a numpy array without coercing its dtype, so
    // the returned event-time array can match the dtype of `y`.
    let y_array: Bound<'py, PyAny> = py.import("numpy")?.getattr("asarray")?.call1((y,))?;

    let (y_dtype, n) = {
        let untyped = y_array.cast::<PyUntypedArray>()?;
        if untyped.ndim() != 1 {
            return Err(PyValueError::new_err(format!(
                "y should be 1-dimensional, got an array with {} dimensions",
                untyped.ndim(),
            )));
        }
        (untyped.dtype(), untyped.len())
    };

    if y_observed.len() != n || weight.as_ref().is_some_and(|w| w.len() != n) {
        return Err(PyValueError::new_err(
            "all arguments should be array-like with the same length",
        ));
    }
    if let Some(ref w) = weight {
        if !w.as_array().iter().all(|&w| w >= 0.0 && w.is_finite()) {
            return Err(PyValueError::new_err(
                "weights should be nonnegative and finite",
            ));
        }
    }

    let y_observed_view = y_observed.as_array();
    let weight_view = weight.as_ref().map(|w| w.as_array());

    fn run<'py, T: TimeValue>(
        py: Python<'py>,
        y_array: Bound<'py, PyAny>,
        y_observed: ArrayView1<'_, bool>,
        weight: Option<ArrayView1<'_, f64>>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyArray1<f64>>)> {
        let typed: Bound<'py, PyArray1<T>> = y_array.cast_into()?;
        let readonly = typed.readonly();
        let jumps = kaplan_meier_jumps::<T>(readonly.as_array(), y_observed, weight);
        let (times, survival): (Vec<T>, Vec<f64>) = jumps.into_iter().unzip();
        Ok((
            PyArray1::from_vec(py, times).into_any(),
            PyArray1::from_vec(py, survival),
        ))
    }

    macro_rules! dispatch {
        ($($T:ty),+ $(,)?) => {$(
            if y_dtype.is_equiv_to(&dtype::<$T>(py)) {
                return run::<$T>(py, y_array, y_observed_view, weight_view);
            }
        )+};
    }

    dispatch!(f64, f32, i64, i32, i16, i8, u64, u32, u16, u8);

    Err(PyValueError::new_err(format!(
        "unsupported dtype for y: {}",
        y_dtype.str()?,
    )))
}

mod park;

#[allow(clippy::type_complexity)]
#[pyfunction(
    signature = (
        x,
        y,
        y_observed,
        centers = None,
        epsilon = 1e-4,
        parallel = false,
    )
)]
fn fit_park<'py>(
    py: Python<'py>,
    x: PyArrayLike1<'py, f64, AllowTypeChange>,
    y: PyArrayLike1<'py, f64, AllowTypeChange>,
    y_observed: PyArrayLike1<'py, bool, AllowTypeChange>,
    centers: Option<PyArrayLike1<'py, f64, AllowTypeChange>>,
    epsilon: f64,
    parallel: bool,
) -> PyResult<(
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    let mut storage = None;
    let covariate = maybe_allocate(&x, &mut storage);
    let mut storage = None;
    let response = maybe_allocate(&y, &mut storage);
    let mut storage = None;
    let observed = maybe_allocate(&y_observed, &mut storage);
    let mut storage = None;
    let mut maybe_centers = None;
    if let Some(centers) = &centers {
        maybe_centers = Some(maybe_allocate(centers, &mut storage));
    }

    let ParkFitResult {
        fit,
        bucket_centers,
        times,
    } = park::fit_park(
        covariate,
        response,
        observed,
        maybe_centers,
        epsilon,
        parallel,
    );

    Ok((
        PyArray1::from_vec(py, fit)
            .reshape((bucket_centers.len(), times.len()))
            .unwrap(),
        PyArray1::from_vec(py, bucket_centers),
        PyArray1::from_vec(py, times),
    ))
}

fn get_slice<'a, D: Dimension>(
    view: &'a ArrayView<f64, D>,
    owned: &'a mut Option<Array<f64, D>>,
) -> &'a [f64] {
    view.as_slice().unwrap_or_else(|| {
        *owned = Some(view.to_owned());
        owned.as_ref().unwrap().as_slice().unwrap()
    })
}

fn maybe_allocate_view<'a, A: Copy, D: Dimension + 'a>(
    array: &'a ArrayView<'a, A, D>,
    storage: &'a mut Option<Vec<A>>,
) -> &'a [A] {
    array.as_slice().unwrap_or_else(|| {
        *storage = Some(array.iter().copied().collect());
        storage.as_deref().unwrap()
    })
}

fn maybe_allocate<'a, T: Copy + Element, D: Dimension>(
    array: &'a PyArrayLike<T, D, AllowTypeChange>,
    storage: &'a mut Option<Vec<T>>,
) -> &'a [T] {
    array.as_slice().unwrap_or_else(|_| {
        *storage = Some(array.as_array().into_iter().copied().collect());
        storage.as_deref().unwrap()
    })
}

/// Isotonic distributional regression (IDR) is a powerful nonparametric technique
/// for the estimation of distributions of a binary or numeric response variable
/// conditional on numeric or ordinal covariates. IDR assumes that there is a
/// monotone relationship between the response variable and the covariates, where
/// the partial order on the covariate space can be specified by the user, and has
/// no tuning parameters. It can be used to generate calibrated probabilistic
/// weather forecasts from ensemble forecasts and observations, and serve as a
/// benchmark in many other prediction problems.
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<IDR>()?;
    m.add_function(wrap_pyfunction!(isotonic_regression, m)?)?;
    m.add_function(wrap_pyfunction!(kaplan_meier, m)?)?;
    m.add_function(wrap_pyfunction!(fit_park, m)?)?;
    Ok(())
}
