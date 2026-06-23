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
