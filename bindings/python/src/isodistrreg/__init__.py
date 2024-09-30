from importlib.metadata import version

__version__ = version("isodistrreg")

from isodistrreg._core import (
    fit_park,
    IDR,
    isotonic_regression,
    kaplan_meier,
)

try:
    import sklearn
    from isodistrreg.estimator import IsotonicDistributionalRegressor
except ImportError:
    pass
