from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt

def fit_park(
    x: npt.ArrayLike,
    time: npt.ArrayLike,
    event: npt.ArrayLike,
    centers: Optional[npt.ArrayLike],
    epsilon: Optional[float],
    parallel: Optional[bool],
):
    pass

class IDR:
    """
    Fit the Isotonic Distributional Regression (IDR) model.

    Parameters
    ----------
    y:
        1D real array-like of shape (n,).
    X:
        Real array-like of shape (n,) or (n, d). Must have d >= 1. The covariate dimension d influences the accepted
        shapes of covariate arguments going forward.
    y_observed:
        Optional 1D boolean array-like of shape (n,). If provided, True marks observed, False marks censored.
    sample_weight:
        Optional 1D real array-like of shape (n,). Must be non-negative if provided.
    X_order:
        Optional sequence of (name, indices) where indices is a sequence of 0-based column indices.
    y_order:
        Optional string; default "sd". May be None, which results also in "sd".
    decreasing:
        If True, fit is monotone decreasing in responses; default False.
    subsamples:
        Optional integer, how many subsamples to create while subagging (default is 1: a single sample with all data)
    subsample_size:
        Optional float (strictly between 0.0 and 1.0) or integer (at least 1, less than the number of observations) used
        for subagging (default is half the data).
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
    def predict(self, X: npt.ArrayLike) -> npt.NDArray[np.float64]: ...

    """
    Predict the conditional mean (or point prediction) for each row in covariates.

    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    Returns: float64 ndarray.
    """

    def cdf(self, X: npt.ArrayLike) -> npt.NDArray[np.float64]: ...

    """
    Predict conditional CDF values over the model’s internal grid for each row.
    
    The thresholds to which the CDF values correspond can be accessed via the "thresholds" field.

    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    Returns: float64 ndarray of shape (..., t).
    """

    def cdf_at(
        self,
        X: npt.ArrayLike,
        y: npt.ArrayLike,
    ) -> npt.NDArray[np.float64]: ...

    """
    Predict the CDF at the specified covariates and thresholds.
    
    X: Array-like of any shape if fit argument had shape (n,), or array-like of shape (..., d) if fit argument had
        shape (n, d).
    y: Array-like of any shape broadcastable with X if fit argument had shape (n,), or broadcastable with all but the 
        last dimension of X if fit argument had shape (n, d).
    Returns: float64 ndarray.
    """

    def cdf_grid(
        self,
        X: npt.ArrayLike,
        y: npt.ArrayLike,
    ) -> npt.NDArray[np.float64]: ...

    """
    Predict the CDF at a grid of **sorted** covariate and **sorted** response values.
    
    X: 1D array-like of shape (m,) if fit argument had shape (n,), or array-like of shape (m, d) if fit argument had
        shape (n, d).
    y: 1D array-like of shape (t,).
    Returns: float64 ndarray of shape (m, t).
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

    def X(self) -> npt.NDArray[np.float64]: ...
    def thresholds(self) -> npt.NDArray[np.float64]: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

def isotonic_regression(
    y: npt.ArrayLike,
    X: Optional[npt.ArrayLike] = None,
    sample_weight: Optional[npt.ArrayLike] = None,
    decreasing: bool = False,
) -> npt.NDArray[np.float64]: ...
def kaplan_meier(
    times: npt.ArrayLike,
    events: npt.ArrayLike,
    weights: Optional[npt.ArrayLike] = None,
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]: ...
