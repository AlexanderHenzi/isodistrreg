# isodistrreg: Python bindings

<p align="center"><strong>Available on <a href="https://pypi.org/project/isodistrreg/">PyPI</a></strong></p>

Python bindings for Isotonic Distributional Regression (IDR) and Survival-IDR
(S-IDR), built with [PyO3](https://pyo3.rs) and
[maturin](https://www.maturin.rs). See the [main
README](https://github.com/AlexanderHenzi/isodistrreg) for background.

## Installation

```
pip install isodistrreg
```

Pre-built wheels are provided for Linux, macOS, and Windows on CPython 3.13+.
On other platforms pip builds from the source distribution, which requires a
[Rust toolchain](https://rustup.rs).

The scikit-learn–compatible estimator is available through an optional extra:

```
pip install "isodistrreg[scikit-learn]"
```

## Examples

We use numpy to set up a toy problem, but plain python lists can also be used.

### Precision

CDF outputs (`cdf`, `cdf_at`, `cdf_grid`) are always `np.float32`. Covariates and
thresholds are stored at the dtype of the input arrays passed to `IDR(...)`:
`float32` inputs stay `float32`, `float64` inputs stay `float64`, and `.X` /
`.thresholds` return zero-copy views at the storage dtype. See the docstring on
`IDR.cdf` and `IDR.from_cdfs` in `_core.pyi` for the full rules.

### Example 1: Covariate 1-dimensional, outcome censored

```python
import numpy as np
# we visualize with matplotlib
from matplotlib import pyplot as plt
from isodistrreg import IDR

# Generate an instance of increasing conditional CDFs with censoring, and a one-dimensional covariate
n = 500
rng = np.random.default_rng(seed=123)
x = rng.uniform(size=n)
y = x + rng.uniform(size=n)
c = x + rng.uniform(size=n)

t = np.minimum(y, c)
d = y <= c

# Fit the IDR / S-IDR model, censoring is indicated by "False"
fit = IDR(t, x, d)

# The same fit from one structured array with a numeric time field and a
# boolean event field (any names, either order, e.g. scikit-survival's layout)
td = np.empty(n, dtype=[("event", "?"), ("time", "f8")])
td["time"], td["event"] = t, d
fit = IDR(td, x)

# Sorted and deduplicated covariates and thresholds are available
sorted_x = fit.X
sorted_y = fit.thresholds

# Estimate and plot the complete distributional estimate
cdf_for_each_x = fit.cdf(sorted_x)
def plot_cdfs_nicely(cdfs, centers, times):
    """Plot a picture with the right scale; plt.imshow is simpler, but spacing along the axis is not to scale"""
    plt.pcolormesh(
        [centers[0] - (centers[1] - centers[0]) / 2]
            + list((centers[1:] + centers[:-1]) / 2)
            + [centers[-1] + (centers[-1] - centers[-2]) / 2],
        list(times) + [times[-1] + (times[-1] - times[0]) / len(times)],
        cdfs.T,
        vmin=0.0,
        vmax=1.0,
    )
plot_cdfs_nicely(cdf_for_each_x, sorted_x, sorted_y)
plt.colorbar()

# Estimate and plot the mean
mean_for_each_x = fit.predict(sorted_x)
# Due to censoring, the estimated sub-CDF may not always have a mean
plt.plot(sorted_x, sorted_x + 0.5, color="lightblue", label="true mean")
plt.plot(sorted_x, mean_for_each_x, color="red", label="mean")

# Estimate and plot quantiles
probabilities = np.array([0.2, 0.8])
quantiles = fit.quantile(sorted_x[:, np.newaxis], probabilities)
plt.plot(sorted_x, quantiles, label=[f"{p} quantile" for p in probabilities])

plt.legend(loc="lower right")
plt.show()
```

### Example 2: Covariate 3-dimensional
```python
import numpy as np
from isodistrreg import IDR

## toy data (3-dimensional covariate)
X = np.column_stack([np.arange(1, 5)] * 3)
y = np.array([1, 0, 2, 2])

## fit
idr_fit = IDR(X = X, y = y)

## get CDF for new x at all relevant thresholds
new_x = np.array([[1, 1, 1], [1.5, 1.5, 1.5]])
idr_fit.cdf(new_x) # (one CDF per row = per x)

## broadcasting
idr_fit.cdf_at(new_x, 0) # (evaluate CDF at 0 for all covariates)
idr_fit.cdf_at(new_x, [0,1]) # (evaluate CDF for x1 at 0, x2 at 1)
idr_fit.cdf_at(new_x, np.column_stack([[1, 2, 3], [0, 1, 2]])) # (CDF at 1,2,3 for x1, and 0,1,2 for x2)

## same for quantiles
idr_fit.quantile(new_x, 0.5)
idr_fit.quantile(new_x, [0.25, 0.5])
idr_fit.quantile(new_x, np.column_stack([[0.25, 0.5, 0.75], [0.1, 0.2, 0.3]]))

# Fast grid evaluation for 1-dimensional covariate
X = np.arange(5)
y = np.arange(5)
idr_fit = IDR(X = X, y = y)
idr_fit.cdf_grid(X, y) # (CDF at all covariate-threshold-combinations)
```

### Example 3: scikit-learn estimator

`IsotonicDistributionalRegressor` wraps `IDR` in the scikit-learn estimator
API, so it works with pipelines, cross-validation and model selection. Model
settings are constructor parameters; per-sample data (`sample_weight`,
`y_observed`) goes to `fit`.

```python
import numpy as np
from sklearn.model_selection import cross_val_score
from isodistrreg import IsotonicDistributionalRegressor

rng = np.random.default_rng(seed=123)
X = rng.uniform(size=(500, 1))
y = X[:, 0] + rng.uniform(size=500)

# All parameters are optional; these are the defaults except random_state.
model = IsotonicDistributionalRegressor(
    covariate_order=None,   # partial order on the columns of X, e.g. [("sd", [0, 1])]
    response_order="sd",    # or "hazard"
    decreasing=False,
    subsamples=None,        # set to subag: fit on random subsamples and average
    random_state=0,         # seeds the subsample draws
).fit(X, y)

# Point predictions (conditional means) and cross-validated R^2 of them
model.predict(X[:3])                      # shape (3,)
cross_val_score(model, X, y, cv=5)

# The distributional estimate: CDF on the fitted threshold grid, at chosen
# thresholds, or its quantiles. Every row is evaluated at every threshold.
model.cdf(X[:3])                          # shape (3, len(model.thresholds_))
model.cdf_at(X[:3], [0.5, 1.0, 1.5])      # shape (3, 3)
model.quantile(X[:3], [0.1, 0.5, 0.9])    # shape (3, 3)
model.quantile(X[:3], 0.5)                # shape (3,)

# Right-censored responses (S-IDR): pass the event indicator to fit, or pass
# y as a structured array with a numeric time and a boolean event field
observed = rng.uniform(size=500) < 0.7
survival = IsotonicDistributionalRegressor().fit(X, y, y_observed=observed)
survival.cdf_at(X[:3], 1.0)               # shape (3,)
```

Under scikit-learn's metadata routing, request the per-sample arguments as
usual, e.g. `model.set_fit_request(sample_weight=True)`.

## References

Henzi, A., Ziegel, J. and Gneiting, T. (2021). Isotonic distributional
regression. *J R Stat Soc Series B*, 83: 963–993.
<https://doi.org/10.1111/rssb.12450>

Henzi, A., Moesching, A. and Duembgen, L. (2022). Accelerating the
Pool-Adjacent-Violators Algorithm for Isotonic Distributional Regression.
*Methodol Comput Appl Probab*.
<https://doi.org/10.1007/s11009-022-09937-2>

Bladt, M., Henzi, A., van den Heuvel, B. and Ziegel, J. (2026). Survival
Isotonic Distributional Regression. *arXiv preprint* arXiv:2608.02914.
<https://arxiv.org/abs/2608.02914>
