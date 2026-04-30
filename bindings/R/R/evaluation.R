#' Cumulative distribution function (CDF) of IDR or raw forecasts
#'
#' @description
#' Evaluate the the cumulative distribution function (CDF) of IDR predictions or
#' of unprocessed forecasts in a \code{data.frame}.
#'
#' @usage
#' cdf(predictions, thresholds)
#'
#' @param predictions either an object of class \code{idr} (output of
#'   \code{\link{predict.idrfit}}), or a \code{data.frame} of numeric variables.
#'   In the latter case, the CDF is computed using the empirical distribution of
#'   the variables in \code{predictions}.
#' @param thresholds numeric vector of thresholds at which the CDF will be
#'   evaluated.
#'
#' @details
#' The CDFs are considered as piecewise constant stepfunctions: If \code{x} are
#' the points where the IDR fitted CDF (or the empirical distribution of the
#' forecasts) has jumps and \code{p} the corresponding CDF values, then for
#' \code{x[i] <= x < x[i + 1]}, the CDF at \code{x} is \code{p[i]}.
#'
#' @return
#' A matrix of probabilities giving the evaluated CDFs at the given thresholds,
#' one column for each threshold.
#'
#' @seealso
#' \code{\link{predict.idrfit}} \code{\link{qpred}}, \code{\link{bscore}}
#'
#' @export
#'
#' @examples
#'
#' data("rain")
#'
#' ## Postprocess HRES forecast using data of 3 years
#'
#' X <- rain[1:(3 * 365), "HRES", drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idr(y = y, X = X)
#'
#' ## Compute probability of precipitation given that the HRES forecast is
#' ## 0 mm, 0.5 mm or 1 mm
#'
#' predictions <- predict(fit, data = data.frame(HRES = c(0, 0.5, 1)))
#' 1 - cdf(predictions, thresholds = 0)
cdf <- function(predictions, thresholds) {
  UseMethod("cdf")
}

#' cdf method for class 'idr'
#'
#' @method cdf idr
#' @rdname cdf
#' @export
cdf.idr <- function(predictions, thresholds) {
  if (!is.vector(thresholds, "numeric")) {
    stop("'thresholds' must be a numeric vector")
  }

  predictions$external_ptr$cdf_at(
    predictions$predict_covariates,
    thresholds
  )
}

#' cdf method for class 'data.frame'
#'
#' @method cdf data.frame
#' @rdname cdf
#' @importFrom stats ecdf
#' @export
cdf.data.frame <- function(predictions, thresholds) {
  if (!is.vector(thresholds, "numeric")) {
    stop("'thresholds' must be a numeric vector")
  }
  if (!all(sapply(predictions, is.numeric))) {
    stop("'predictions' contains non-numeric variables")
  }
  predictions <- asplit(data.matrix(predictions), 1)
  cdf0 <- function(data) {
    stats::ecdf(data)(thresholds)
  }
  cdfVals <- lapply(predictions, cdf0)
  do.call(rbind, unname(cdfVals))
}


#' Quantile function of IDR or raw forecasts
#'
#' @description
#' Evaluate the the quantile function of IDR predictions or of unprocessed
#' forecasts in a \code{data.frame}.
#'
#' @usage
#' qpred(predictions, quantiles)
#'
#' @param predictions either an object of class \code{idr} (output of
#'   \code{\link{predict.idrfit}}), or a \code{data.frame} of numeric variables.
#'   In the latter case, quantiles are computed using the empirical distribution
#'   of the variables in \code{predictions}.
#' @param quantiles numeric vector of desired quantiles.
#'
#' @details
#' The quantiles are defined as lower quantiles, that is,
#' \deqn{
#'   q(u) = inf(x: cdf(x) >= u),
#' }
#' except for \deqn{u = 0} when the lower endpoint of the support is returned.
#'
#' @return
#' A matrix of forecasts for the desired quantiles, one column per quantile.
#'
#' @seealso
#' \code{\link{predict.idrfit}}, \code{\link{cdf}}, \code{\link{qscore}}
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## Postprocess HRES forecast using data of 3 years
#'
#' X <- rain[1:(3 * 365), "HRES", drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idr(y = y, X = X)
#'
#' ## Compute 95%-quantile forecast given that the HRES forecast is
#' ## 2.5 mm, 5 mm or 10 mm
#'
#' predictions <- predict(fit, data = data.frame(HRES = c(2.5, 5, 10)))
#' qpred(predictions, quantiles = 0.95)
qpred <- function(predictions, quantiles) {
  UseMethod("qpred")
}

#' qpred method for class 'idr'
#'
#' @method qpred idr
#'
#' @importFrom stats stepfun
#' @importFrom utils tail
#' @rdname qpred
#' @export
qpred.idr <- function(predictions, quantiles) {
  # Check input
  if (!is.vector(quantiles, "numeric") || min(quantiles) < 0 ||
        max(quantiles) > 1) {
    stop("quantiles must be a numeric vector with entries in [0,1]")
  }

  predictions$external_ptr$quantile(
    predictions$predict_covariates,
    quantiles
  )
}


#' qpred method for class 'data.frame'
#'
#' @method qpred data.frame
#'
#' @rdname qpred
#' @importFrom stats quantile
#' @export
qpred.data.frame <- function(predictions, quantiles) {
  # Check input
  if (!is.vector(quantiles, "numeric") || min(quantiles) < 0 ||
        max(quantiles) > 1) {
    stop("quantiles must be a numeric vector with entries in [0,1]")
  }
  if (!all(sapply(predictions, is.numeric))) {
    stop("'predictions' contains non-numeric variables")
  }

  predictions <- asplit(data.matrix(predictions), 1)
  q0 <- function(x) {
    stats::quantile(x, probs = quantiles, type = 1)
  }
  qVals <- lapply(predictions, q0)
  unname(do.call(rbind, qVals))
}

#' Quantile scores for IDR or raw forecasts
#'
#' @description
#' Computes quantile scores of IDR quantile predictions or of quantile
#' predictions from raw forecasts in a \code{data.frame}.
#'
#' @usage
#' qscore(predictions, quantiles, y)
#'
#' @inheritParams qpred
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions, or of length 1. In the latter case, \code{y} will be used
#'   for all predictions.
#'
#' @details
#' The quantile score of a forecast \emph{x} for the \emph{u}-quantile is
#' defined as
#' \deqn{
#' 2(1{x > y} - u)(x - y),
#' }
#' where \emph{y} is the observation. For \emph{u = 1/2}, this equals the mean
#' absolute error of the median forecast.
#'
#' @return
#' A matrix of the quantile scores for the desired quantiles, one column per
#' quantile.
#'
#' @seealso
#' \code{\link{predict.idrfit}}, \code{\link{qpred}}
#'
#' @references
#' Gneiting, T. and Raftery, A. E. (2007), 'Strictly proper scoring rules,
#' prediction, and estimation', Journal of the American Statistical Association
#' 102(477), 359-378
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## Postprocess HRES forecast using data of 3 years
#'
#' X <- rain[1:(3 * 365), "HRES", drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idr(y = y, X = X)
#'
#' ## Compute mean absolute error of the median postprocessed forecast using
#' ## data of the next 2 years (out-of-sample predictions) and compare to raw
#' ## HRES forecast
#'
#' data <- rain[(3 * 365 + 1):(5 * 365), "HRES", drop = FALSE]
#' obs <- rain[(3 * 365 + 1):(5 * 365), "obs"]
#'
#' predictions <- predict(fit, data = data)
#' idrMAE <- mean(qscore(predictions, 0.5, obs))
#' rawMAE <- mean(qscore(data, 0.5, obs))
#'
#' c("idr" = idrMAE, "raw" = rawMAE)
qscore <- function(predictions, quantiles, y) {
  if (!is.vector(y, "numeric")) {
    stop("y must be a numeric vector")
  }
  if (!is.vector(quantiles)) {
    quantiles <- c(quantiles)
  }
  predicted <- qpred(predictions, quantiles)
  if (!length(y) %in% c(1, nrow(predicted))) {
    stop("y must have length 1 or the same length as the predictions")
  }
  qsVals <- sweep(
    x = predicted,
    MARGIN = 1,
    STATS = y
  )
  2 * qsVals * sweep(qsVals > 0, MARGIN = 2, STATS = quantiles)
}

#' Brier score for forecast probability of threshold exceedance
#'
#' @description Computes the Brier score of forecast probabilities for exceeding
#' given thresholds.
#'
#' @usage bscore(predictions, thresholds, y)
#'
#' @inheritParams cdf
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions, or of length 1. In the latter case, \code{y} will be used for
#'   all predictions.
#'
#' @details The Brier score for the event of exceeding a given threshold
#' \emph{z} is defined as \deqn{ (1\{y > z\} - P(y > z))^2 } where \emph{y} is
#' the observation and \emph{P(y > z)} the forecast probability for exceeding
#' the threshold \code{z}.
#'
#' @return A matrix of the Brier scores for the desired thresholds, one column
#' per threshold.
#'
#' @seealso \code{\link{predict.idrfit}}, \code{\link{cdf}}
#'
#' @references
#' Gneiting, T. and Raftery, A. E. (2007), 'Strictly proper scoring rules,
#' prediction, and estimation', Journal of the American Statistical Association
#' 102(477), 359-378
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## Postprocess HRES forecast using data of 3 years
#'
#' X <- rain[1:(3 * 365), "HRES", drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idr(y = y, X = X)
#'
#' ## Compute Brier score for postprocessed probability of precipitation
#' ## forecast using data of the next 2 years (out-of-sample predictions)
#'
#' data <- rain[(3 * 365 + 1):(5 * 365), "HRES", drop = FALSE]
#' obs <- rain[(3 * 365 + 1):(5 * 365), "obs"]
#' predictions <- predict(fit, data = data)
#' score <- bscore(predictions, thresholds = 0, y = obs)
#'
#' mean(score)
bscore <- function(predictions, thresholds, y) {
  # Check input
  if (!is.vector(y, "numeric")) {
    stop("obs must be a numeric vector")
  }
  predicted <- cdf(predictions, thresholds)
  if (!(lY <- length(y)) %in% c(1, nrow(predicted))) {
    stop("y must have length 1 or the same length as the predictions")
  }

  if (lY == 1) {
    bsVals <- sweep(
      x = predicted,
      MARGIN = 2,
      STATS = (y <= thresholds)
    )
  } else {
    bsVals <- predicted - outer(y, thresholds, "<=")
  }
  bsVals^2
}

#' Continuous ranked probability score (CRPS)
#'
#' @description Computes the CRPS of IDR or raw forecasts.
#'
#' @usage crps(predictions, y)
#'
#' @param predictions either an object of class \code{idr} (output of
#'   \code{\link{predict.idrfit}}), or a \code{data.frame} of numeric variables.
#'   In the latter case, the CRPS is computed using the empirical distribution
#'   of the variables in \code{predictions}.
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions, or of length 1. In the latter case, \code{y} will be used for
#'   all predictions.
#'
#' @details
#' This function uses adapted code taken from the function \code{crps_edf} of
#' the \pkg{scoringRules} package.
#'
#' @return A vector of CRPS values.
#'
#' @seealso \code{\link{predict.idrfit}}
#'
#' @references
#' Jordan A., Krueger F., Lerch S. (2018). "Evaluating Probabilistic
#' Forecasts with scoringRules." Journal of Statistical Software. Forthcoming.
#'
#' Gneiting, T. and Raftery, A. E. (2007), 'Strictly proper scoring rules,
#' prediction, and estimation', Journal of the American Statistical Association
#' 102(477), 359-378
#'
#' @export
#'
#' @examples
#' data("rain")
#'
#' ## Postprocess HRES forecast using data of 3 years
#'
#' X <- rain[1:(3 * 365), "HRES", drop = FALSE]
#' y <- rain[1:(3 * 365), "obs"]
#'
#' fit <- idr(y = y, X = X)
#'
#' ## Compute CRPS of postprocessed HRES forecast using data of the next 2 years
#' ## (out-of-sample predictions)
#'
#' data <- rain[(3 * 365 + 1):(5 * 365), "HRES", drop = FALSE]
#' obs <- rain[(3 * 365 + 1):(5 * 365), "obs"]
#' predictions <- predict(fit, data = data)
#' idrCrps <- crps(predictions, y = obs)
#'
#' ## Compare this to CRPS of the raw ensemble of all forecasts (high
#' ## resolution, control and 50 perturbed ensemble forecasts)
#'
#' rawData <- rain[(3 * 365 + 1):(5 * 365), c("HRES", "CTR", paste0("P", 1:50))]
#' rawCrps <- crps(rawData, y = obs)
#'
#' c("idr_HRES" = mean(idrCrps), "raw_all" = mean(rawCrps))
crps <- function(predictions, y) {
  UseMethod("crps")
}

#' crps method for class 'irdpred'
#'
#' @method crps idr
#' @rdname crps
#' @export
crps.idr <- function(predictions, y) {
  # Check input
  if (!is.vector(y, "numeric")) {
    stop("obs must be a numeric vector")
  }
  if (length(y) != 1 && length(y) != nrow(predictions$cdf)) {
    stop("y must have length 1 or the same length as the predictions")
  }

  p <- predictions$points
  cdfs <- predictions$cdf
  w <- cbind(cdfs[, 1], t(apply(cdfs, 1, diff)))

  vapply(
    seq_along(y),
    function(i) {
      2 * sum(w[i, ] * ((y[i] < p) - cdfs[i, ] + 0.5 * w[i, ]) * (p - y[i]))
    },
    numeric(1)
  )
}

#' crps method for class 'data.frame'
#'
#' @method crps data.frame
#' @rdname crps
#' @export
crps.data.frame <- function(predictions, y) {
  # Check input
  if (!is.vector(y, "numeric")) {
    stop("obs must be a numeric vector")
  }
  n <- nrow(predictions)
  if (length(y) != 1 && length(y) != n) {
    stop("y must have length 1 or the same length as the predictions")
  }

  x <- lapply(split(data.matrix(predictions), seq_len(n)), sort)
  m <- ncol(predictions)
  a <- seq.int(0.5 / m, 1 - 0.5 / m, length.out = m)
  crps0 <- function(y, x, a) {
    2 / m * sum(((y < x) - a) * (x - y))
  }
  mapply(crps0,
    y = y,
    x = x,
    MoreArgs = list(a = a)
  )
}


#' Probability integral transform (PIT)
#'
#' @description Computes the probability integral transform (PIT) of IDR or raw
#'   forecasts.
#'
#' @usage pit(predictions, y, randomize = TRUE, seed = NULL)
#'
#' @param predictions either an object of class \code{idr} (output of
#'   \code{\link{predict.idrfit}}), or a \code{data.frame} of numeric variables.
#'   In the latter case, the PIT is computed using the empirical distribution of
#'   the variables in \code{predictions}.
#' @param y a numeric vector of obervations of the same length as the number of
#'   predictions.
#' @param randomize PIT values should be randomized at discontinuity points of
#'   the predictive CDF (e.g. at zero for precipitation forecasts). Set \code{
#'   randomize = TRUE} to randomize.
#' @param seed argument to \code{set.seed} for random number generation (if
#'   \code{randomize} is \code{TRUE}).
#'
#' @return Vector of PIT values.
#'
#' @seealso \code{\link{predict.idrfit}}
#'
#' @references
#'
#' Gneiting, T., Balabdaoui, F. and Raftery, A. E. (2007), 'Probabilistic
#' forecasts, calibration and sharpness', Journal of the Royal Statistical
#' Society: Series B (Statistical Methodology) 69(2), 243-268.
#'
#' @export
#'
#' @importFrom stats runif
#'
#' @examples
#' data("rain")
#' require("graphics")
#'
#' ## Postprocess HRES forecast using data of 4 years
#'
#' X <- rain[1:(4 * 365), "HRES", drop = FALSE]
#' y <- rain[1:(4 * 365), "obs"]
#'
#' fit <- idr(y = y, X = X)
#'
#' ## Assess calibration of the postprocessed HRES forecast using data of next 4
#' ## years and compare to calibration of the raw ensemble
#'
#' data <- rain[(4 * 365 + 1):(8 * 365), "HRES", drop = FALSE]
#' obs <- rain[(4 * 365 + 1):(8 * 365), "obs"]
#' predictions <- predict(fit, data = data)
#' idrPit <- pit(predictions, obs, seed = 123)
#'
#' rawData <- rain[(4 * 365 + 1):(8 * 365), c("HRES", "CTR", paste0("P", 1:50))]
#' rawPit <- pit(rawData, obs, seed = 123)
#'
#' hist(idrPit,
#'   xlab = "Probability Integral Transform",
#'   ylab = "Density", freq = FALSE, main = "Postprocessed HRES"
#' )
#' hist(rawPit,
#'   xlab = "Probability Integral Transform",
#'   ylab = "Density", freq = FALSE, main = "Raw ensemble"
#' )
pit <- function(predictions,
                y,
                randomize = TRUE,
                seed = NULL) {
  UseMethod("pit")
}


#' pit method for class 'idr'
#'
#' @method pit idr
#' @rdname pit
#' @importFrom stats stepfun
#' @export
pit.idr <- function(predictions,
                    y,
                    randomize = TRUE,
                    seed = NULL) {
  # Check input
  if (!is.vector(y, "numeric")) {
    stop("'y' must be a numeric vector")
  }
  if (length(y) != nrow(predictions$cdf)) {
    stop("'y' must have the same length as the predictions")
  }

  fitted_thresholds <- predictions$points
  pit0 <- function(cdf, y) {
    # Evaluated CDF (stepfun)
    stats::stepfun(x = fitted_thresholds, y = c(0, cdf))(y)
  }
  cdfs <- predictions$cdf
  pitVals <- vapply(
    seq_len(nrow(cdfs)),
    function(i) pit0(cdfs[i, ], y[i]),
    numeric(1)
  )

  if (randomize) {
    # Randomization: Find small epsilon to compute left limit of CDF
    eps <- if (length(fitted_thresholds) > 1) {
      min(diff(fitted_thresholds))
    } else {
      1
    }

    # Left limits F(y-): evaluate at y - eps/2, row by row
    lowerPitVals <- vapply(
      seq_len(nrow(cdfs)),
      function(i) pit0(cdfs[i, ], y[i] - eps / 2),
      numeric(1)
    )

    if (!is.null(seed)) {
      set.seed(seed)
    }

    idx <- which(lowerPitVals < pitVals)
    pitVals[idx] <- stats::runif(length(idx),
      min = lowerPitVals[idx],
      max = pitVals[idx]
    )
  }

  pitVals
}

#' pit method for class 'data.frame'
#'
#' @method pit data.frame
#' @rdname pit
#' @importFrom stats ecdf
#' @export
pit.data.frame <- function(predictions,
                           y,
                           randomize = TRUE,
                           seed = NULL) {
  # Check input
  if (!is.vector(y, "numeric")) {
    stop("'y' must be a numeric vector")
  }
  n <- nrow(predictions)
  if (length(y) != n) {
    stop("'y' must have the same length as the predictions")
  }
  if (!all(sapply(predictions, is.numeric))) {
    stop("'predictions' contains non-numeric variables")
  }

  pit0 <- function(data, y) {
    stats::ecdf(data)(y)
  }
  pred <- unname(split(data.matrix(predictions), seq_len(n)))
  # Update split(data.matrix(...)) to asplit
  pitVals <- mapply(pit0, data = pred, y = y)
  if (randomize) {
    # Randomization: Find small epsilon to compute left limit of CDF
    eps <- apply(predictions, 1, stats::dist, method = "manhattan")
    eps <- min(c(eps[eps > 0], 1))
    lowerPitVals <- mapply(pit0, data = pred, y = y - eps / 2)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    sel <- lowerPitVals < pitVals
    if (any(sel)) {
      pitVals[sel] <- stats::runif(
        sum(sel),
        min = lowerPitVals[sel],
        max = pitVals[sel]
      )
    }
  }
  pitVals
}

#' Plot IDR predictions
#'
#' @description Plot an IDR predictive CDF.
#'
#' @method plot idr
#'
#' @param x object of class \code{idr} (output of
#'   \code{\link{predict.idrfit}}).
#' @param index index of the prediction in \code{x} for which a plot is desired.
#' @param bounds whether the bounds should be plotted or not (see
#'   \code{\link{predict.idrfit}} for details about the meaning of the bounds).
#' @param col.cdf color of the predictive CDF.
#' @param col.bounds color of the bounds.
#' @param lty.cdf linetype of the predictive CDF.
#' @param lty.bounds linetype of the CDF bounds.
#' @param main main title.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param ... further arguments to \code{\link{plot.stepfun}} or
#'   \code{\link{plot}}.
#'
#' @return
#' The data based on which the plot is drawn (returned invisible).
#'
#' @seealso \code{\link{predict.idrfit}}
#'
#' @export
#' @importFrom stats stepfun
#' @importFrom graphics plot
#'
#' @examples
#' data("rain")
#' require("graphics")
#'
#' ## Postprocess HRES and CTR forecast using data of 2 years
#'
#' X <- rain[1:(2 * 365), c("HRES", "CTR"), drop = FALSE]
#' y <- rain[1:(2 * 365), "obs"]
#'
#' ## Fit IDR and plot the predictive CDF when the HRES forecast is 1 mm and
#' ## CTR is 0 mm
#'
#' fit <- idr(y = y, X = X)
#' pred <- predict(fit, data = data.frame(HRES = 1, CTR = 0))
#' plot(pred)
plot.idr <- function(x,
                     index = 1,
                     bounds = TRUE,
                     col.cdf = "black",
                     col.bounds = "blue",
                     lty.cdf = 1,
                     lty.bounds = 3,
                     xlab = "Threshold",
                     ylab = "CDF",
                     main = "IDR predictive CDF",
                     ...) {
  cdf <- x$cdf[index, ]

  graphics::plot(
    stepfun(x = x$points, y = c(0, cdf)),
    xlab = xlab,
    ylab = ylab,
    do.points = FALSE,
    col = col.cdf,
    lty = lty.cdf,
    main = main,
    ...
  )

  graphics::abline(
    h = c(0, 1),
    col = "gray",
    lty = 3
  )

  invisible(cdf)
}
