#' Compute IDR predictions with (su)bagging
#'
#' @description Computes IDR predictions with bootstrap aggregating (bagging)
#' or subsample aggregation (subagging).
#' 
#' @usage idrbag(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)), orders =
#'   c("comp" = 1), stoch = "sd", pars = osqpSettings(verbose = FALSE, eps_abs =
#'   1e-5, eps_rel = 1e-5, max_iter = 10000L), progress = TRUE, newdata, 
#'   digits = 3, interpolation = "linear", asplitAvail = TRUE, b, p, 
#'   replace = FALSE, grid = NULL)
#' 
#' @param newdata \code{data.frame} containing variables with which to
#'   predict. Ordered factor variables are converted to numeric for computation,
#'   so ensure that the factor levels are identical in \code{newdata} and in
#'   \code{X}.
#' @param digits number of decimal places for the predictive CDF.
#' @param interpolation interpolation method for univariate data. Default is 
#'   \code{"linear"}. Any other argument will select midpoint interpolation (see 
#'   'Details' in \code{\link{predict.idrfit}}). Has no effect for multivariate
#'   IDR.
#' @param asplitAvail use \code{\link[base]{asplit}} for splitting arrays
#'   (default is \code{TRUE}). Set to \code{FALSE} for R Versions < 3.6, where
#'   \code{asplit} is not available.
#' @param b number of (su)bagging samples.
#' @param p size of (su)bagging samples relative to training data.
#' @param replace draw samples with (\code{TRUE}) or without (\code{FALSE}) 
#'     replacement?
#' @param grid grid on which the predictive CDFs are evaluated. Default are
#'     the unique values of \code{y}.
#' @inheritParams idr
#' 
#' @details 
#' This function draws \code{b} times a random subsample of size
#' \code{ceiling(nrow(X)*p)}) from the training data, fits IDR to each
#' subsample, computes predictions for the new data supplied in \code{newdata},
#' and averages the predictions derived from the \code{b} subsamples. There are
#' no default values for \code{b} and \code{p}.
#' 
#' @return
#' A list of predictions, see \code{\link{predict.idrfit}}).
#' 
#' @export
idrbag <- function(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)),
  orders = c("comp" = 1), stoch = "sd", pars = osqpSettings(verbose = FALSE,
  eps_abs = 1e-5, eps_rel = 1e-5, max_iter = 10000L), progress = TRUE, newdata,
  digits = 3, interpolation = "linear", asplitAvail = TRUE, b, p, 
  replace = FALSE, grid = NULL) {
    
  if (!is.vector(y, mode = "numeric")) 
    stop("'y' must be a numeric vector")
  N <- length(y)
  if (!is.numeric(b) | length(b) != 1 | !(as.integer(b) == b) | b < 1)
    stop("'b' must be a positive integer smaller than length(y)")
  if (!is.numeric(p) | length(p) != 1 | p <= 0 | p >= 1)
    stop("'p' must be a number in (0,1)")
  if (!isTRUE(replace) & !isFALSE(replace))
    stop("'replace' must be TRUE or FALSE")
  if (!isTRUE(progress) & !isFALSE(progress))
    stop("'replace' must be TRUE or FALSE")
  if (is.null(grid)) {
    grid <- sort(unique(y))
  } else {
    if (!is.vector(grid, "numeric"))
      stop("'grid' must be a numeric vector or NULL")
    grid <- sort(grid)
  }
  if (!is.data.frame(newdata))
    stop("'newdata' must be a data.frame")
  
  n <- ceiling(p * N)
  m <- length(grid)
  preds <- matrix(nrow = nrow(newdata), ncol = m, 0)
  if (progress) {
    pb <- utils::txtProgressBar(max = b)
    for (i in seq_len(b)) {
      utils::setTxtProgressBar(pb, i)
      s <- sample(N, n, replace = replace)
      fit <- idr(y = y[s], X = X[s, , drop = FALSE], groups = groups,
        orders = orders, stoch = stoch, pars = pars, progress = FALSE)
      preds <- preds + cdf(predict(object = fit, data = newdata, digits = digits,
        interpolation = interpolation, asplitAvail = asplitAvail), grid)
    }
    close(pb)
  } else {
    for (i in seq_len(b)) {
      s <- sample(N, n, replace = replace)
      fit <- idr(y = y[s], X = X[s, , drop = FALSE], groups = groups,
        orders = orders, stoch = stoch, pars = pars, progress = FALSE)
      preds <- preds + cdf(predict(object = fit, data = newdata, digits = digits,
        interpolation = interpolation, asplitAvail = asplitAvail), grid)
    }
  }
  preds <- asplit(round(preds / b, digits), 1)
  preds <- lapply(
    X = preds,
    FUN = function(dat) {
      sel <- c(dat[1] > 0, dat[-1] > dat[-m])
      data.frame(points = grid[sel], cdf = dat[sel])
    }
  )
  structure(preds, class = "idr", incomparables = integer(0))
}