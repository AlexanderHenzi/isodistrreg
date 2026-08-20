# isodistrreg 0.6.0

* Multivariate (partial order) fits now use a quadratic-program solver built
  into the package instead of OSQP. Installation no longer needs CMake or a C
  toolchain, which removes a class of installation failure on platforms where
  CMake is absent or misconfigured. `pars` keeps the same options with the same
  meanings and the same defaults, so existing code fits the same models.
* Fixed two failures on hazard-rate-order fits of larger datasets: one aborted
  the fit outright, the other silently reported full convergence. Both are now
  reported through `diagnostic$convergence_fraction`.
* `diagnostic$convergence_fraction` now reports what it documents. It could
  previously only ever be `1`, so a fit that did not converge was
  indistinguishable from one that did.
* Rust users: `Config`'s `osqp_settings` field is now `solver_settings`, and its
  type is owned by the crate rather than by the solver. Python users: the
  `"osqp_settings"` settings key is now `"solver_settings"`; the keys inside it
  are unchanged.

# isodistrreg 0.5.2

* Faster fitting and prediction, especially for censored and larger datasets.
* Fixed an incorrect fit when data mixed negative and positive zero (`-0`/`0`).
* Fixed a prediction rounding bug that could break quantile monotonicity in the
  covariate.
* Updated Rust and R build dependencies.

# isodistrreg 0.5.1

* `idrbag()` now returns a proper `"idr"` prediction object and validates its
  `newdata`/`grid` inputs.
* Fixed `crps()` for CDFs with one or two jump points and for scalar `y`, and
  `pit()` to randomize only at actual jump points.
* Fixed ordered factor covariates and by-name column matching in `idr()`,
  `idrbag()` and `predict.idrfit()`, which now also returns in-sample
  predictions when `data` is omitted.
* `pars` is now merged with the documented defaults, and
  `isotonic_regression()` validates its inputs.
* Invalid inputs now raise targeted R errors instead of Rust panics, and the
  progress bar is written to stderr.
* Removed the inoperative `bounds`, `col.bounds` and `lty.bounds` arguments of
  `plot.idr()`, validated `index`, and corrected the documentation.

# isodistrreg 0.5.0

* switched core computations to f32 for performance

# isodistrreg 0.4.2

* improve performance for the partial order case under censoring

# isodistrreg 0.4.1

* improve performance for the total order case under censoring

# isodistrreg 0.4.0

* add support for right-censored outcomes in `idr` and `idrbag` for partial
  orders

# isodistrreg 0.3.0

* add support for right-censored outcomes in `idr` and `idrbag` for total
  orders

# isodistrreg 0.2.0

* add weights to `idr` function

# isodistrreg 0.1.0

* isodistrreg is now on CRAN: https://CRAN.R-project.org/package=isodistrreg
* removed argument `asplitAvail` from the functions. The package now requires
  R version 3.6 or newer
