from importlib.metadata import version
from typing import TYPE_CHECKING

__version__ = version("isodistrreg")

from isodistrreg._core import (
    fit_park,
    IDR,
    isotonic_regression,
    kaplan_meier,
)

if TYPE_CHECKING:
    from isodistrreg.estimator import IsotonicDistributionalRegressor

# The scikit-learn estimator is left out of ``__all__`` on purpose: it needs the
# optional scikit-learn dependency, so ``from isodistrreg import *`` must not
# try to resolve it.
__all__ = ["IDR", "fit_park", "isotonic_regression", "kaplan_meier"]


def __getattr__(name: str):
    # Resolve the estimator lazily: importing it pulls in scikit-learn, and a
    # missing scikit-learn should be reported by name at the point of use rather
    # than silently hiding the class.
    if name == "IsotonicDistributionalRegressor":
        try:
            from isodistrreg.estimator import IsotonicDistributionalRegressor
        except ModuleNotFoundError as e:
            missing = (e.name or "").partition(".")[0]
            if missing not in ("sklearn", "joblib"):
                raise
            raise ImportError(
                "isodistrreg.IsotonicDistributionalRegressor requires scikit-learn; "
                'install it with `pip install "isodistrreg[scikit-learn]"`'
            ) from e
        return IsotonicDistributionalRegressor
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
