# isodistrreg: Isotonic Distributional Regression

Rust, Python and R implementations of **Isotonic Distributional Regression
(IDR)** and its **Survival-IDR (S-IDR)** extension for right-censored data.

## What is IDR?

IDR is a nonparametric method for estimating the *conditional distribution* of
a numeric or binary response given numeric or ordinal covariates. Instead of
predicting a single point (like a mean or a quantile), it returns the full
conditional CDF, from which means, quantiles, exceedance probabilities, and
proper scoring rules can all be derived.

Compared to other distributional regression techniques, IDR has a few
properties that make it appealing as a baseline or as a calibration layer:

- **No tuning parameters.** The fit is determined entirely by the data and the
  user-specified partial order on the covariate space.
- **Monotone by construction.** IDR assumes a monotone relationship between
  covariates and response, so the estimated CDFs respect stochastic ordering.
- **Calibrated and sharp.** The estimator is the unique solution to a
  constrained optimization that is optimal under a wide class of proper scoring
  rules (CRPS, Brier score, quantile loss, ...).

A motivating use case is post-processing ensemble weather forecasts: IDR turns
a raw ensemble into a calibrated probabilistic forecast without distributional
assumptions, and serves as a strong benchmark in many other prediction tasks.

### S-IDR: extension to right-censored outcomes

In survival analysis, the outcome is often only observed up to a censoring
time. **S-IDR** extends IDR to this setting: given `(t, d)` pairs where `t` is
the observed time and `d` indicates whether `t` is the true event time or a
right-censoring time, S-IDR estimates the conditional sub-distribution of the
event time. All the usual IDR queries (CDF, quantiles, mean where defined)
remain available.

## Language bindings

The same Rust core powers all three interfaces. Pick the one that matches your
stack:

| Language   | Path                                            |
|------------|-------------------------------------------------|
| **Python** | [`bindings/python`](bindings/python/README.md)  |
| **R**      | [`bindings/R`](bindings/R/README.md)            |
| **Rust**   | [`isodistrreg`](isodistrreg/README.md)          |

Each subdirectory's README contains installation instructions and worked
examples.

## References

Henzi, A., Ziegel, J.F. and Gneiting, T. (2021). Isotonic distributional
regression. *J R Stat Soc Series B*, 83: 963–993.
<https://doi.org/10.1111/rssb.12450>

Henzi, A., Moesching, A. and Duembgen, L. (2022). Accelerating the
Pool-Adjacent-Violators Algorithm for Isotonic Distributional Regression.
*Methodol Comput Appl Probab*.
<https://doi.org/10.1007/s11009-022-09937-2>
