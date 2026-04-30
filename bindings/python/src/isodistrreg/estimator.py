import logging
from typing import Self, cast

import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted, validate_data

import numpy.typing as npt

from isodistrreg import IDR


class IsotonicDistributionalRegressor(RegressorMixin, BaseEstimator):
    def fit(self, X: npt.ArrayLike, y, **kwargs) -> Self:
        X, y = validate_data(self, X, y)
        X = cast(np.ndarray, X)
        self.n_features_in_ = X.shape[-1]

        # Parse optional arguments
        y_observed = kwargs.pop("y_observed", None)
        if y_observed is not None:
            y_observed = check_array(y_observed)
            if not np.isin(y_observed, [0, 1]).all():
                raise ValueError(
                    "y_observed must contain only 0/1 or True/False values"
                )
        sample_weight = kwargs.pop("sample_weight", None)
        if sample_weight is not None:
            sample_weight = check_array(sample_weight)
        if self.n_features_in_ > 1:
            covariate_order = kwargs.pop("covariate_order", None)
            if covariate_order is not None:
                seen = [False] * self.n_features_in_
                for kind, ids in covariate_order:
                    if kind not in ("comp", "sd", "icx"):
                        raise ValueError(
                            f"covariate_order kind must be 'comp', 'sd', or 'icx', got '{kind}'"
                        )
                    for i in ids:
                        if not isinstance(i, int) or i >= self.n_features_in_:
                            raise ValueError(
                                f"covariate_order index must be an int < {self.n_features_in_}, got {i!r}"
                            )
                        if seen[i]:
                            raise ValueError(f"duplicate column '{i}'")
                        seen[i] = True
                del seen
        else:
            covariate_order = None
        response_order = kwargs.pop("response_order", None)
        if response_order is not None and response_order not in ("sd", "hazard"):
            raise ValueError(
                f"response_order must be 'sd' or 'hazard', got '{response_order}'"
            )
        decreasing = kwargs.pop("decreasing", False)
        if decreasing not in (True, False):
            raise ValueError(f"decreasing must be True or False, got {decreasing!r}")
        subsamples = kwargs.pop("subsamples", None)
        subsample_size = kwargs.pop("subsample_size", None)

        if kwargs:
            # Unused arguments
            logging.info(f"Unused kwargs {list(kwargs.keys())}")

        self._fit = IDR(
            y,
            X,
            y_observed,
            sample_weight,
            covariate_order,
            response_order,
            decreasing,
            subsamples,
            subsample_size,
        )
        self.X_ = self._fit.X
        self.thresholds_ = self._fit.thresholds

        return self

    def predict(self, X: npt.ArrayLike, **kwargs) -> np.typing.NDArray:
        check_is_fitted(self)
        targets = validate_data(self, X, reset=False)

        if kwargs:
            # Unused arguments
            logging.info(f"Unused kwargs {list(kwargs.keys())}")

        return self._fit.predict(targets)

    def cdf(self, X: npt.ArrayLike, **kwargs) -> np.typing.NDArray:
        check_is_fitted(self)
        targets = validate_data(self, X, reset=False)

        if kwargs:
            # Unused arguments
            logging.info(f"Unused kwargs {list(kwargs.keys())}")

        return self._fit.cdf(targets)

    def cdf_at(
        self, data: npt.ArrayLike, thresholds: npt.ArrayLike, **kwargs
    ) -> np.typing.NDArray:
        check_is_fitted(self)
        # Check thresholds as y values
        targets = validate_data(self, data, thresholds, reset=False)

        if kwargs:
            # Unused arguments
            logging.info(f"Unused kwargs {list(kwargs.keys())}")

        return self._fit.cdf_at(targets, thresholds)

    def quantiles(
        self, data: npt.ArrayLike, quantiles: npt.ArrayLike, **kwargs
    ) -> np.typing.NDArray:
        check_is_fitted(self)
        targets = validate_data(self, data, reset=False)
        if self.n_features_in_ != 1:
            raise ValueError(
                "quantiles are not defined for multidimensional covariates X"
            )

        quantiles = np.asarray(quantiles, dtype=np.float64)
        if not np.all((0.0 <= quantiles) & (quantiles <= 1.0)):
            raise ValueError("quantiles must be between 0 and 1")

        # Predict the lower quantile by default
        upper = kwargs.pop("upper", False)

        if kwargs:
            # Unused arguments
            logging.info(f"Unused kwargs {list(kwargs.keys())}")

        return self._fit.quantile(targets, quantiles, upper)
