"""scikit-learn estimator wrapping :class:`isodistrreg.IDR`."""

import operator
from numbers import Integral, Real
from typing import Any, Self, cast

import numpy as np
import numpy.typing as npt
from joblib import effective_n_jobs
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils import check_random_state
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.validation import (
    _check_sample_weight,
    check_array,
    check_consistent_length,
    check_is_fitted,
    column_or_1d,
    validate_data,
)

from isodistrreg import IDR
from isodistrreg._core import split_censored_outcome


class IsotonicDistributionalRegressor(RegressorMixin, BaseEstimator):
    """Isotonic distributional regression (IDR) as a scikit-learn estimator.

    IDR estimates the full conditional distribution of a response ``y`` given
    covariates ``X`` under the assumption that the distribution is monotone in
    the covariates: larger covariates give stochastically larger responses (or
    smaller ones with ``decreasing=True``). No tuning parameters are needed for
    the fit itself; the parameters below select the ordering assumptions and
    optional subsample aggregation. Right-censored responses are supported
    through ``y_observed`` at fit time, or by passing ``y`` as a structured
    array holding both times and event indicators (survival IDR, S-IDR).

    Parameters
    ----------
    covariate_order : list of (str, list of int), default=None
        Partial order on the covariate space. Each entry is a pair
        ``(kind, column_indices)``, where ``kind`` is ``"comp"`` (componentwise
        order), ``"sd"`` (stochastic dominance of the empirical distribution
        of the listed columns) or ``"icx"`` (increasing convex order). Columns
        not listed in any entry are ordered componentwise. ``None`` orders all
        columns componentwise. Supplying an order selects the partial-order
        solver even for a single covariate.
    response_order : {"sd", "hazard"}, default="sd"
        Stochastic order imposed on the conditional distributions along the
        covariate order: usual stochastic dominance, or hazard rate order.
    decreasing : bool, default=False
        Fit distributions that are stochastically decreasing in the covariates
        instead of increasing.
    subsamples : int, default=None
        Number of subsamples to fit and average (subagging). ``None`` fits the
        model once on all data. Any value, including ``1``, draws random
        subsamples and therefore uses ``random_state``.
    subsample_size : int or float, default=None
        Size of each subsample: an absolute count if an int, a fraction of the
        training set if a float in ``(0, 1]``. Defaults to half the training
        set without replacement and to the full training set with replacement.
        Only meaningful together with ``subsamples``.
    replace : bool, default=False
        Draw subsamples with replacement (bootstrap) instead of without
        (subagging).
    settings : dict, default=None
        Solver settings forwarded to :class:`isodistrreg.IDR`. Currently the
        only key is ``"solver_settings"``, used by the partial-order solver.
    random_state : int, RandomState instance or None, default=None
        Controls the subsample draws when ``subsamples`` is set. Pass an int
        for reproducible fits. Ignored without subsampling.
    n_jobs : int, default=None
        Number of threads used to fit the subsamples in parallel. ``None``
        means 1, ``-1`` means all processors. Ignored without subsampling.

    Attributes
    ----------
    X_unique_ : ndarray of shape (n_unique, n_features)
        Sorted, deduplicated training covariates; the covariate grid on which
        the fitted conditional CDFs are stored. Stored at the dtype of the
        training covariates.
    thresholds_ : ndarray of shape (n_thresholds,)
        Sorted, deduplicated training responses; the grid of thresholds on
        which the fitted conditional CDFs are stored. Stored at the dtype of
        the training responses.
    n_features_in_ : int
        Number of features seen during :term:`fit`.
    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    isodistrreg.IDR : The underlying model, with a lower-level array interface.

    Notes
    -----
    Conditional CDF values are computed in single precision (``float32``)
    regardless of the input dtypes; see the ``isodistrreg`` README.

    The inherited :meth:`score` is the coefficient of determination of
    :meth:`predict`, i.e. of the estimated conditional means. It ignores the
    rest of the estimated distribution and is only a crude summary of a
    distributional fit.

    References
    ----------
    Henzi, A., Ziegel, J. and Gneiting, T. (2021). Isotonic distributional
    regression. *J R Stat Soc Series B*, 83: 963-993.

    Bladt, M., Henzi, A., van den Heuvel, B. and Ziegel, J. (2026). Survival
    Isotonic Distributional Regression. arXiv:2608.02914.

    Examples
    --------
    >>> import numpy as np
    >>> from isodistrreg import IsotonicDistributionalRegressor
    >>> X = np.array([[1.0], [2.0], [3.0], [4.0]])
    >>> y = np.array([10.0, 20.0, 15.0, 30.0])
    >>> idr = IsotonicDistributionalRegressor().fit(X, y)
    >>> idr.predict([[2.5]])
    array([17.5])
    >>> idr.cdf_at([[2.5]], [15.0, 20.0])
    array([[0.5, 1. ]], dtype=float32)
    >>> idr.quantile([[2.5]], [0.25, 0.75])
    array([[15., 20.]])
    """

    _parameter_constraints: dict[str, list[Any]] = {
        "covariate_order": [list, tuple, None],
        "response_order": [StrOptions({"sd", "hazard"})],
        "decreasing": ["boolean"],
        "subsamples": [Interval(Integral, 1, None, closed="left"), None],
        "subsample_size": [
            Interval(Integral, 1, None, closed="left"),
            Interval(Real, 0, 1, closed="right"),
            None,
        ],
        "replace": ["boolean"],
        "settings": [dict, None],
        "random_state": ["random_state"],
        "n_jobs": [Integral, None],
    }

    def __init__(
        self,
        *,
        covariate_order: list[tuple[str, list[int]]] | None = None,
        response_order: str = "sd",
        decreasing: bool = False,
        subsamples: int | None = None,
        subsample_size: int | float | None = None,
        replace: bool = False,
        settings: dict[str, Any] | None = None,
        random_state: int | np.random.RandomState | None = None,
        n_jobs: int | None = None,
    ) -> None:
        self.covariate_order = covariate_order
        self.response_order = response_order
        self.decreasing = decreasing
        self.subsamples = subsamples
        self.subsample_size = subsample_size
        self.replace = replace
        self.settings = settings
        self.random_state = random_state
        self.n_jobs = n_jobs

    def fit(
        self,
        X: npt.ArrayLike,
        y: npt.ArrayLike,
        sample_weight: npt.ArrayLike | None = None,
        y_observed: npt.ArrayLike | None = None,
    ) -> Self:
        """Fit the model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training covariates.
        y : array-like of shape (n_samples,)
            Responses. With censoring, the observed times; alternatively a
            structured array with a boolean event-indicator field and a
            numeric time field, in either order and under any names, such as
            ``[("time", "f8"), ("event", "?")]`` or scikit-survival's
            ``[("event", "?"), ("time", "f8")]``. ``y_observed`` is then
            omitted.
        sample_weight : array-like of shape (n_samples,), default=None
            Non-negative observation weights. ``None`` weights all
            observations equally. Weights are processed in single precision.
        y_observed : array-like of shape (n_samples,), default=None
            Right-censoring indicator: ``True`` (or ``1``) marks an observed
            response, ``False`` (or ``0``) a censored one. ``None`` treats all
            responses as observed unless ``y`` is a structured array carrying
            its own indicator.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        self._validate_params()
        # A structured y carries the event indicator; unpack it before
        # scikit-learn's validation, which only accepts a numeric y.
        y, y_observed = split_censored_outcome(y, y_observed)
        X, y = validate_data(self, X, y, y_numeric=True)

        if sample_weight is not None:
            sample_weight = _check_sample_weight(
                sample_weight, X, ensure_non_negative=True
            )
        observed = None
        if y_observed is not None:
            # Shape and length are checked here for scikit-learn's messages;
            # the values (boolean, or numeric 0/1 only) are checked by the core.
            observed = column_or_1d(
                check_array(y_observed, ensure_2d=False, dtype=None), warn=True
            )
            check_consistent_length(X, observed)

        # The core draws its subsamples from a u64 seed; derive one from the
        # scikit-learn random_state only when subsampling is in effect, so a
        # plain fit never consumes random numbers.
        if self.subsamples is None:
            seed = None
        else:
            random_state = check_random_state(self.random_state)
            seed = int(random_state.randint(np.iinfo(np.int32).max))

        self._fit = IDR(
            y,
            X,
            y_observed=observed,
            sample_weight=sample_weight,
            X_order=self._covariate_order_entries(),
            y_order=self.response_order,
            decreasing=self.decreasing,
            subsamples=self.subsamples,
            subsample_size=self.subsample_size,
            replace=self.replace,
            settings=self.settings,
            seed=seed,
            n_jobs=effective_n_jobs(self.n_jobs),
        )
        self.X_unique_ = self._fit.X
        self.thresholds_ = self._fit.thresholds

        return self

    def _covariate_order_entries(self) -> list[tuple[str, list[int]]] | None:
        """Normalise ``covariate_order`` to what the core accepts.

        Only the structure is checked here; the order kinds, index bounds and
        duplicate columns are validated by the core.
        """
        if self.covariate_order is None:
            return None
        entries = []
        for entry in self.covariate_order:
            try:
                kind, columns = entry
            except (TypeError, ValueError):
                raise ValueError(
                    "each covariate_order entry must be a (kind, column_indices) "
                    f"pair, got {entry!r}"
                ) from None
            if not isinstance(kind, str):
                raise ValueError(f"covariate_order kind must be a str, got {kind!r}")
            indices: list[int] = []
            for column in columns:
                if isinstance(column, bool) or not isinstance(column, Integral):
                    raise ValueError(
                        f"covariate_order column indices must be ints, got {column!r}"
                    )
                indices.append(operator.index(column))
            entries.append((kind, indices))
        return entries

    def predict(self, X: npt.ArrayLike) -> npt.NDArray[np.floating]:
        """Predict the conditional mean for each row of ``X``.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Covariates.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,)
            Estimated conditional means, at the dtype of the training
            responses.
        """
        covariates = self._check_covariates(X)
        return self._fit.predict(covariates)

    def cdf(self, X: npt.ArrayLike) -> npt.NDArray[np.float32]:
        """Estimate the conditional CDF of each row of ``X`` on the threshold grid.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Covariates.

        Returns
        -------
        cdf : ndarray of shape (n_samples, n_thresholds)
            ``cdf[i, j]`` is the estimated probability that the response is at
            most ``thresholds_[j]`` given covariates ``X[i]``.
        """
        covariates = self._check_covariates(X)
        return self._fit.cdf(covariates)

    def cdf_at(
        self, X: npt.ArrayLike, thresholds: npt.ArrayLike
    ) -> npt.NDArray[np.float32]:
        """Estimate the conditional CDF of each row of ``X`` at given thresholds.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Covariates.
        thresholds : float or array-like of shape (n_thresholds,)
            Thresholds at which to evaluate every row's conditional CDF.

        Returns
        -------
        cdf : ndarray of shape (n_samples,) or (n_samples, n_thresholds)
            ``cdf[i, j]`` is the estimated probability that the response is at
            most ``thresholds[j]`` given covariates ``X[i]``. The trailing
            axis is dropped for a scalar ``thresholds``.
        """
        covariates = self._check_covariates(X)
        thresholds = self._check_grid(thresholds, "thresholds")
        # Insert a broadcast axis so each row is evaluated at every threshold.
        out = self._fit.cdf_at(covariates[:, None, :], np.atleast_1d(thresholds))
        return out[:, 0] if thresholds.ndim == 0 else out

    def quantile(
        self, X: npt.ArrayLike, q: npt.ArrayLike, *, upper: bool = False
    ) -> npt.NDArray[np.floating]:
        """Estimate conditional quantiles for each row of ``X``.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Covariates.
        q : float or array-like of shape (n_quantiles,)
            Probability levels in ``[0, 1]`` at which to evaluate every row's
            conditional quantile function.
        upper : bool, default=False
            Return the upper quantile (the largest threshold at which the CDF
            equals the level) instead of the lower one where the two differ.

        Returns
        -------
        quantiles : ndarray of shape (n_samples,) or (n_samples, n_quantiles)
            ``quantiles[i, j]`` is the estimated ``q[j]``-quantile of the
            response given covariates ``X[i]``. The trailing axis is dropped
            for a scalar ``q``.
        """
        covariates = self._check_covariates(X)
        q = self._check_grid(q, "q")
        if not np.all((0.0 <= q) & (q <= 1.0)):
            raise ValueError("q must be between 0 and 1")
        # Insert a broadcast axis so each row is evaluated at every level.
        out = self._fit.quantile(covariates[:, None, :], np.atleast_1d(q), upper)
        return out[:, 0] if q.ndim == 0 else out

    def _check_covariates(self, X: npt.ArrayLike) -> npt.NDArray[Any]:
        """Validate prediction-time covariates against the fitted model."""
        check_is_fitted(self)
        return cast(npt.NDArray[Any], validate_data(self, X, reset=False))

    @staticmethod
    def _check_grid(values: npt.ArrayLike, name: str) -> npt.NDArray[np.float64]:
        values = np.asarray(values, dtype=np.float64)
        if values.ndim > 1:
            raise ValueError(f"{name} must be a scalar or 1-D array-like")
        return values
