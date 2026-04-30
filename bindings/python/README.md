# isodistrreg: Python bindings

Python bindings for Isotonic Distributional Regression (IDR) and Survival-IDR
(S-IDR), built with [PyO3](https://pyo3.rs) and
[maturin](https://www.maturin.rs). See the [main
README](https://github.com/AlexanderHenzi/isodistrreg) for background and
references.

## Installation

```bash
pip install isodistrreg
```

## Parallelism

The package supports optional parallelism via [rayon](https://docs.rs/rayon).
To install with parallelism enabled, build from source with the `parallel`
Cargo feature:

```bash
maturin develop --features parallel
```

When built with parallelism, you can control the number of threads at runtime
using the `RAYON_NUM_THREADS` environment variable:

```bash
# Use 4 threads
RAYON_NUM_THREADS=4 python my_script.py

# Use a single thread (disables parallelism at runtime)
RAYON_NUM_THREADS=1 python my_script.py
```

If `RAYON_NUM_THREADS` is not set, rayon defaults to using all available CPU
cores.

## Examples

We use numpy to set up a toy problem, but plain python lists can also be used.

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

# Sorted and deduplicated covariates and thresholds are available
sorted_x = fit.covariates
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
quantiles = fit.quantile(sorted_x[:, np.newaxis], probabilities[np.newaxis])
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
