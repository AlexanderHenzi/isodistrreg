#' Compute IDR predictions with (su)bagging
#'
#' @description Computes IDR predictions with bootstrap aggregating (bagging)
#' or subsample aggregation (subagging).
#'
#' @usage idrbag(y, X, y_observed = NULL, weights = NULL, decreasing = FALSE,
#'   groups = setNames(rep(1, ncol(X)), colnames(X)), orders = c("comp" = 1),
#'   stoch = "sd", pars = list(verbose = FALSE, eps_abs = 1e-5,
#'   eps_rel = 1e-5, max_iter = 10000L), n_jobs = 1, progress = TRUE, newdata,
#'   digits = NULL, interpolation = "linear", b, p, replace = FALSE,
#'   grid = NULL, seed = NULL)
#'
#' @param newdata \code{data.frame} containing variables with which to
#'   predict. Ordered factor variables are converted to numeric for computation,
#'   so ensure that the factor levels are identical in \code{newdata} and in
#'   \code{X}.
#' @param digits removed functionality, parameter kept for backwards
#'   compatibility but ignored with warning: number of decimal places for the
#'   predictive CDF, useful to keep the solution small across covariates.
#' @param interpolation interpolation method for univariate data. Default is
#'   \code{"linear"}. Any other argument will select midpoint interpolation (see
#'   'Details' in \code{\link{predict.idrfit}}). Has no effect for multivariate
#'   IDR.
#' @param b number of (su)bagging samples.
#' @param p size of (su)bagging samples relative to training data.
#' @param replace draw samples with (\code{TRUE}, \code{1}) or without
#'     (\code{FALSE}, \code{0}) replacement?
#' @param grid grid on which the predictive CDFs are evaluated. Default are
#'     the unique values of \code{y}.
#' @param seed integer seed for the random number generator. Only relevant
#'   when (su)bagging is active.
#' @param n_jobs number of worker threads used to fit the individual subsamples
#'   in parallel. Only relevant when (su)bagging is active. Default is \code{1}
#'   (serial execution).
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
#' A list of predictions, see \code{\link{predict.idrfit}}.
#'
#' @export
idrbag <- function(y,
                   X,
                   y_observed = NULL,
                   weights = NULL,
                   decreasing = FALSE,
                   groups = setNames(rep(1, ncol(X)), colnames(X)),
                   orders = c("comp" = 1),
                   stoch = "sd",
                   pars = list(
                     verbose = FALSE,
                     eps_abs = 1e-5,
                     eps_rel = 1e-5,
                     max_iter = 10000L
                   ),
                   n_jobs = 1,
                   progress = TRUE,
                   newdata = NULL,
                   digits = NULL,
                   interpolation = "linear",
                   b,
                   p,
                   replace = FALSE,
                   grid = NULL,
                   seed = NULL) {
  inputs <- validate(
    y,
    X,
    y_observed,
    weights,
    decreasing,
    groups,
    orders,
    stoch,
    pars,
    progress,
    seed
  )
  if (!is.numeric(b) || length(b) != 1 || !(as.integer(b) == b) || b < 1) {
    stop("'b' must be a positive integer smaller than length(y)")
  }
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("'p' must be a number in (0,1)")
  }
  if (!is.numeric(n_jobs) || length(n_jobs) != 1 ||
        !(as.integer(n_jobs) == n_jobs) || n_jobs < 1) {
    stop("'n_jobs' must be a positive integer")
  }
  n_jobs <- as.integer(n_jobs)
  if (isTRUE(replace == 1)) replace <- TRUE
  if (isTRUE(replace == 0)) replace <- FALSE
  if (isTRUE(progress == 1)) progress <- TRUE
  if (isTRUE(progress == 0)) progress <- FALSE
  if (!isTRUE(replace) && !isFALSE(replace)) {
    stop("'replace' must be TRUE/FALSE or 1/0")
  }
  if (!isTRUE(progress) && !isFALSE(progress)) {
    stop("'progress' must be TRUE/FALSE or 1/0")
  }
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data.frame")
    }
    if (ncol(newdata) != ncol(X)) {
      stop("'newdata' must have the same columns as 'X'")
    }
  }

  # Pack the column indices into the orders list
  order_group <- list()
  for (i in seq_along(inputs$orders)) {
    order_kind <- names(inputs$orders)[i]
    group_id <- inputs$orders[[i]]
    members <- which(inputs$groups == group_id)
    indices <- match(names(members), names(inputs$X))
    order_group <- append(order_group, setNames(list(indices), order_kind))
  }

  # Call the Rust function
  external_ptr <- IDR$fit(
    y              = inputs$y,
    X              = inputs$X_formatted,
    y_observed     = inputs$y_observed,
    sample_weight  = inputs$weights,
    x_order        = order_group,
    y_order        = inputs$stoch,
    decreasing     = inputs$decreasing,
    subsamples     = b,
    subsample_size = p,
    replace        = replace,
    settings       = inputs$pars,
    seed           = inputs$seed,
    n_jobs         = n_jobs,
    show_progress  = inputs$progress
  )

  if (!is.null(newdata)) {
    X_eval <- t(newdata)
  } else {
    X_eval <- inputs$X_formatted
  }
  if (!is.null(grid)) {
    cdf <- external_ptr$cdf_at(X_eval, grid)
    points <- grid
  } else {
    cdf <- external_ptr$cdf(X_eval)
    points <- external_ptr$thresholds()
  }

  structure(
    list(
      y = inputs$y,
      X = inputs$X,
      # takes care of duplicate obs as opposed to simply taking the backing cdfs
      cdf = cdf,
      points = points,
      weights = inputs$weights,
      response_unique = external_ptr$thresholds(),
      groups = inputs$groups,
      orders = inputs$orders,
      diagnostic = external_ptr$diagnostic(),
      external_ptr = external_ptr
    ),
    class = "idrfit"
  )
}
