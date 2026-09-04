from typing import Any, Dict, Optional, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt

def fit_park(
    x: npt.ArrayLike,
    time: npt.ArrayLike,
    event: Optional[npt.ArrayLike] = None,
    centers: Optional[npt.ArrayLike] = None,
    epsilon: float = 1e-4,
    parallel: bool = False,
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    ``time`` may be a structured array holding both the times and the event
    indicators, as described for ``IDR``; ``event`` is then omitted.
    """

class IDR:
    """
    Fit the Isotonic Distributional Regression (IDR) model.

    Parameters
    ----------
    y:
        1D real array-like of shape (n,). A right-censored outcome may also be passed as a structured array
        with two fields, a boolean event indicator and a numeric time, in either order and under any names —
        ``[("time", "f8"), ("event", "?")]`` or scikit-survival's ``[("event", "?"), ("time", "f8")]``. The
        fields are told apart by dtype; ``y_observed`` must then be omitted.
    X:
        Real array-like of shape (n,) or (n, d). Must have d >= 1. The covariate dimension d influences the accepted
        shapes of covariate arguments going forward.
    y_observed:
        Optional 1D array-like of shape (n,). If provided, True marks observed, False marks censored. Must be a
        boolean array or a numeric array containing only 0 and 1; any other value is rejected rather than read
        as an event. Omit it when ``y`` is a structured array carrying its own indicator.
    sample_weight:
        Optional 1D real array-like of shape (n,). Must be non-negative if provided. Weights are processed in
        single precision; it is up to the caller to avoid extreme imbalance (as a rule of thumb, no weight
        below ~1e-7 of the total weight).
    X_order:
        Optional sequence of (name, indices) where indices is a sequence of 0-based column indices.
        Note: supplying X_order selects the partial-order solver even for 1-D covariates.
    y_order:
        Optional string, "sd" or "hazard"; default "sd". May be None, which results also in "sd".
    decreasing:
        If True, fit is monotone decreasing in responses; default False.
    subsamples:
        Optional integer, how many subsamples to create while subagging (default is 1: a single sample with all data)
    subsample_size:
        Optional float in (0, 1] or integer (at least 1, at most the number of observations) used
        for subagging. The default is half the data when ``replace=False`` (subagging) and the
        full sample size when ``replace=True`` (the classic bootstrap). The full size without
        replacement would reproduce the plain fit and is rejected when ``subsamples > 1``.
    settings:
        Optional dictionary with settings for the fit, depends on the dimension of X.
    """

    def __init__(
        self,
        y: npt.ArrayLike,
        X: npt.ArrayLike,
        y_observed: npt.ArrayLike | None = ...,
        sample_weight: npt.ArrayLike | None = ...,
        X_order: Sequence[tuple[str, Sequence[int]]] | None = ...,
        y_order: str | None = "sd",
        decreasing: bool = False,
        subsamples: int | None = None,
        subsample_size: float | int | None = None,
        replace: bool = False,
        settings: Dict[str, Any] | None = None,
        seed: int | None = None,
        n_jobs: int = 1,
        progress: bool = False,
    ) -> None: ...
    def predict(
        self, X: npt.ArrayLike
    ) -> Union[npt.NDArray[np.float64], npt.NDArray[np.float32]]: ...

    """
    Predict the conditional mean (or point prediction) for each row in covariates.

    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    Returns: ndarray whose dtype matches the model's response dtype (float64 if y was f64, float32 if y was f32).
    """

    def cdf(self, X: npt.ArrayLike) -> npt.NDArray[np.float32]: ...

    """
    Predict conditional CDF values over the model’s internal grid for each row.

    The thresholds to which the CDF values correspond can be accessed via the "thresholds" field.

    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    Returns: float32 ndarray of shape (..., t).
    """

    def cdf_at(
        self,
        X: npt.ArrayLike,
        y: npt.ArrayLike,
    ) -> npt.NDArray[np.float32]: ...

    """
    Predict the CDF at the specified covariates and thresholds.

    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    y: Array-like of any shape broadcastable with X if fit argument had shape (n,), or broadcastable with all but the
        last dimension of X if fit argument had shape (n, d).
    Returns: float32 ndarray.
    """

    def cdf_grid(
        self,
        X: npt.ArrayLike,
        y: npt.ArrayLike,
    ) -> npt.NDArray[np.float32]: ...

    """
    Predict the CDF at a grid of **sorted** covariate and **sorted** response values
    (sortedness is validated).

    Only supported for models fitted on univariate covariates; multivariate models
    raise ValueError.

    X: 1D array-like of shape (m,), sorted in non-decreasing order.
    y: 1D array-like of shape (t,), sorted in non-decreasing order.
    Returns: float32 ndarray of shape (m, t).
    """

    def quantile(
        self,
        X: npt.ArrayLike,
        q: npt.ArrayLike,
        upper: bool = False,
    ) -> npt.NDArray[np.float64]: ...

    """
    Predict conditional quantiles.

    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    q: Array-like of any shape broadcastable with X if fit argument had shape (n,), or broadcastable with all but the
        last dimension of X if fit argument had shape (n, d).
    Returns: float64 ndarray.
    """

    @property
    def X(self) -> Union[npt.NDArray[np.float64], npt.NDArray[np.float32]]: ...
    @property
    def thresholds(
        self,
    ) -> Union[npt.NDArray[np.float64], npt.NDArray[np.float32]]: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    @classmethod
    def from_cdfs(
        cls,
        cdfs: npt.ArrayLike,
        X: npt.ArrayLike,
        y: npt.ArrayLike,
        global_cdf: npt.ArrayLike | None = ...,
    ) -> "IDR": ...

    """
    Create an IDR model from pre-computed CDFs, covariates, and thresholds.

    cdfs: 2D array-like of shape (n_covariates, n_thresholds). Interpreted as ``np.float32``
        internally — the algorithm always stores CDFs in f32, regardless of how the fit was
        produced. Passing an ``np.float64`` array is allowed but is silently narrowed to
        ``np.float32`` with no warning. The narrowed values are bit-identical to what a
        fresh ``IDR(...)`` fit would store (so there is no precision loss relative to a
        normal fit), but precision is lost relative to the caller's f64 input. Prefer
        ``np.float32`` to make the storage precision explicit at the call site.
    X: 1D or 2D array-like of unique covariates in strictly increasing lexicographic order.
        The storage dtype tracks the input dtype (``np.float32`` or ``np.float64``).
    y: 1D array-like of unique thresholds in strictly increasing order. The storage dtype
        tracks the input dtype (``np.float32`` or ``np.float64``).
    global_cdf: Optional 1D array-like used for incomparable covariates in the
        multivariate case. Same dtype convention as ``cdfs``: ``np.float32`` is preferred,
        ``np.float64`` is silently narrowed.
    Returns: A fitted IDR model.
    """

def isotonic_regression(
    y: npt.ArrayLike,
    X: Optional[npt.ArrayLike] = None,
    sample_weight: Optional[npt.ArrayLike] = None,
    decreasing: bool = False,
    constraints: Optional[npt.ArrayLike] = None,
) -> npt.NDArray[np.float64]: ...
def kaplan_meier(
    y: npt.ArrayLike,
    y_observed: Optional[npt.ArrayLike] = None,
    weight: Optional[npt.ArrayLike] = None,
) -> Tuple[npt.NDArray[Any], npt.NDArray[np.float64]]: ...

"""
Compute the (weighted) Kaplan-Meier estimator of a right-censored sample.

y: Event or censoring times, or a structured array holding both the times and the event indicators as
    described for ``IDR``; ``y_observed`` is then omitted and otherwise required.
Returns: the sorted distinct event times (same dtype as ``y``) and the survival probabilities just after each.
"""

def split_censored_outcome(
    y: npt.ArrayLike,
    y_observed: Optional[npt.ArrayLike] = None,
) -> Tuple[npt.ArrayLike, Optional[npt.ArrayLike]]: ...

"""
Resolve a right-censored outcome into its time and event-indicator arrays.

For a structured ``y`` with a numeric time field and a boolean event field (either order, any names), returns
views of the two fields and requires ``y_observed`` to be omitted. Any other ``y`` is returned unchanged
together with ``y_observed``. Raises ValueError for a structured ``y`` of any other layout.
"""
