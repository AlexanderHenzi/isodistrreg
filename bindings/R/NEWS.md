# isodistrreg (development version)

* `idrbag()` now returns a prediction object of class `"idr"` (usable with
  `cdf()`, `qpred()`, `crps()`, `pit()` and `plot()`), as documented. It also
  validates `newdata` and `grid`, matches `newdata` columns by name, and warns
  when the ignored `digits`/`interpolation` arguments are supplied.
* Fixed `crps()` for fits whose CDFs have only one or two jump points, and for
  scalar observations `y` (previously only the first prediction was scored).
* Fixed `pit()` to randomize only at actual jump points of the predictive CDF;
  values strictly between jumps are now deterministic.
* Fixed ordered factor covariates in `idr()`, `idrbag()` and
  `predict.idrfit()`: they are now converted by factor level, as documented
  (previously the conversion produced `NA`s and an internal error).
* `predict.idrfit()` returns in-sample predictions when `data` is omitted, as
  documented, and matches `data` columns to the training data by name instead
  of by position.
* `pars` options are merged with the documented defaults, so partial lists
  like `pars = list(verbose = TRUE)` work (previously an internal error).
* Invalid inputs (non-finite values, `NA` weights or thresholds, all-zero
  weights, fully censored responses, malformed `groups`/`orders`, duplicated
  column names) now produce targeted R errors instead of internal Rust panics
  or silently wrong results.
* `isotonic_regression()` validates its inputs and accepts integer vectors.
* The progress bar is written to stderr, as documented (previously stdout).
* Removed the inoperative `bounds`, `col.bounds` and `lty.bounds` arguments of
  `plot.idr()`; `index` is now validated.
* The documentation of the `idr()` return value, `dindexm()` and
  `predict.dindexfit()` was corrected to match the actual behavior.

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
