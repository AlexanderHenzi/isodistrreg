#' Fit IDR to training data
#'
#' @description Fits isotonic distributional regression (IDR) to a training
#'   dataset.
#'
#' @usage idr(y, X, y_observed = NULL, weights = NULL, decreasing = FALSE,
#'   groups = setNames(rep(1, ncol(X)), colnames(X)), orders = c("comp" = 1),
#'   stoch = "sd", pars = list(verbose = FALSE, eps_abs = 1e-5,
#'   eps_rel = 1e-5, max_iter = 10000L), progress = TRUE)
#'
#' @param y numeric vector (the response variable).
#' @param X data frame of numeric or ordered factor variables (the regression
#'   covariates).
#' @param y_observed vector of indicators (TRUE or 1 for observed, FALSE or 0
#'   for right-censored). At least one observation must be uncensored. Default
#'   is all observed (\code{rep(TRUE, length(y))}).
#' @param weights vector of finite, non-negative weights (same length as y),
#'   at least one of which must be positive; observations with zero weight are
#'   dropped from the fit. Default is all weights equal to one. Weights are
#'   processed in single precision; it is up to the caller to avoid extreme
#'   imbalance (as a rule of thumb, no weight below ~1e-7 of the total weight).
#' @param decreasing boolean indicating whether \code{y} decreases with \code{X}
#'   (by default, it increases with \code{X}).
#' @param groups named vector of length \code{ncol(X)} denoting groups of
#'   variables that are to be ordered with the same order (see 'Details'). Only
#'   relevant if \code{X} contains more than one variable. The same names as in
#'   \code{X} should be used.
#' @param orders named vector giving for each group in \code{groups} the order
#'   that will be applied to this group. Only relevant if \code{X} contains more
#'   than one variable. The names of \code{orders} give the order, the entries
#'   give the group labels. Available options: \code{"comp"} for componentwise
#'   order, \code{"sd"} for stochastic dominance, \code{"icx"} for increasing
#'   convex order (see 'Details). Default is \code{"comp"} for all variables.
#'   The \code{"sd"} and \code{"icx"} orders can only be used with numeric
#'   variables, but not with ordered factors.
#' @param stoch stochastic order constraint used for estimation. Default is
#'   \code{"sd"} for first order stochastic dominance. Use \code{"hazard"} for
#'   hazard rate order (experimental).
#' @param pars parameters for quadratic programming optimization (only relevant
#'   if \code{X} has more than one column), a list with options "verbose" T / F
#'   (verbosity of solver), "eps_abs" positive float, "eps_rel" positive float,
#'   "max_iter" positive integer.
#' @param progress display a progress bar while fitting (\code{TRUE},
#'   \code{FALSE} or \code{1}, \code{0}). Default is \code{TRUE}; the bar is
#'   written to stderr and is best viewed in an interactive R session.
#'
#' @details This function computes the isotonic distributional regression (IDR)
#'   of a response \emph{y} on on one or more covariates \emph{X}. IDR estimates
#'   the cumulative distribution function (CDF) of \emph{y} conditional on
#'   \emph{X} by monotone regression, assuming that \emph{y} is more likely to
#'   take higher values, as \emph{X} increases. Formally, IDR assumes that the
#'   conditional CDF \eqn{F_{y | X = x}(z)} at each fixed threshold \emph{z}
#'   decreases, as \emph{x} increases, or equivalently, that the exceedance
#'   probabilities for any threshold \code{z} \eqn{P(y > z | X = x)} increase
#'   with \emph{x}.
#'
#'   The conditional CDFs are estimated at each threshold in \code{unique(y)}.
#'   This is the set where the CDFs may have jumps. If \code{X} contains more
#'   than one variable, the CDFs are estimated by solving
#'   \code{length(unique(y))} quadratic programs using osqp. This might take a
#'   while if the training dataset is large.
#'
#'   Use the argument \code{groups} to group \emph{exchangeable} covariates.
#'   Exchangeable covariates are indistinguishable except from the order in
#'   which they are labelled (e.g. ensemble weather forecasts, repeated
#'   measurements under the same measurement conditions).
#'
#'   The following orders are available to perform the monotone regression in
#'   IDR: \itemize{ \item Componentwise order (\code{"comp"}): A covariate
#'   vector \code{x1} is greater than \code{x2} if \code{x1[i] >= x2[i]} holds
#'   for all components \code{i}. This is the \emph{standard order used in
#'   multivariate monotone regression} and \emph{should not be used for
#'   exchangeable variables (e.g. perturbed ensemble forecasts)}. \item
#'   Stochastic dominance (\code{"sd"}): \code{x1} is greater than \code{x2} in
#'   the stochastic order, if the (empirical) distribution of the elements of
#'   \code{x1} is greater than the distribution of the elements of \code{x2} (in
#'   first order) stochastic dominance. The \code{"sd"} order is invariant under
#'   permutations of the grouped variables and therefore \emph{suitable for
#'   exchangeable covariables}. \item Increasing convex order (\code{"icx"}):
#'   The \code{"icx"} order can be used for groups of exchangeable variables. It
#'   should be used if the variables have increasing variability, when their
#'   mean increases (e.g. precipitation forecasts or other variables with
#'   right-skewed distributions). More precisely, \code{"icx"} uses the
#'   increasing convex stochastic order on the empirical distributions of the
#'   grouped variables. }
#'
#' @return An object of class \code{"idrfit"} containing the following
#'   components:
#'
#'   \item{\code{X}}{the training covariates as provided, one row per
#'   observation (in input order, including duplicated rows).}
#'
#'   \item{\code{y}}{numeric vector of the training responses.}
#'
#'   \item{\code{cdf}}{matrix containing the estimated CDFs, one CDF per row,
#'   evaluated at \code{response_unique} (see next point). The CDF in the
#'   \code{i}th row corresponds to the estimated conditional distribution of the
#'   response given the covariates values in \code{X[i,]}.}
#'
#'   \item{\code{weights}}{the observation weights as provided (\code{NULL} if
#'   none were given).}
#'
#'   \item{\code{response_unique}}{the thresholds at which the CDFs in
#'   \code{cdf} are evaluated. The entries in \code{cdf[,j]} are the conditional
#'   CDFs evaluated at \code{response_unique[j]}.}
#'
#'   \item{\code{groups}, \code{orders}}{ the groups and orders used for
#'   estimation.}
#'
#'   \item{\code{diagnostic}}{diagnostics of the CDF estimation. For univariate
#'   fits (total order) this is \code{list(epsilon = )}, a bound on the
#'   precision of the CDF estimation (the maximal downwards-step in the CDF
#'   that has been detected). For multivariate fits (partial order) this is
#'   \code{list(precision = , convergence_fraction = )}, where
#'   \code{convergence_fraction} is the fraction of CDF estimations that
#'   converged before hitting the iteration limit \code{max_iter}. Decrease the
#'   parameters \code{eps_abs} and/or \code{eps_rel} or increase
#'   \code{max_iter} in \code{pars} to improve the precision.}
#'
#'
#' @note The function \code{idr} is only intended for fitting IDR model for a
#'   training dataset and storing the results for further processing, but not
#'   for prediction or evaluation, which is done using the output of
#'   \code{\link{predict.idrfit}}.
#'
#'   The fitted object contains an external pointer to memory managed by the
#'   internal Rust library. It is only valid within the R session that created
#'   it: fits saved with \code{saveRDS} cannot be restored in a new session.
#'
#' @seealso The S3 method \code{\link{predict.idrfit}} for predictions based on
#'   an IDR fit.
#'
#' @export
#' @importFrom stats setNames
#' @importFrom utils modifyList
#'
#' @references Henzi, A., Moesching, A. & Duembgen, L. Accelerating the
#'   Pool-Adjacent-Violators Algorithm for Isotonic Distributional Regression.
#'   Methodol Comput Appl Probab (2022).
#'   https://doi.org/10.1007/s11009-022-09937-2
#'
#' Bladt, M., Henzi, A., van den Heuvel, B. and Ziegel, J. (2026). Survival
#' Isotonic Distributional Regression. arXiv:2608.02914.
#' https://doi.org/10.48550/arXiv.2608.02914
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020).
#' OSQP: An operator splitting solver for quadratic programs. Mathematical
#' Programming Computation, 1-36.
#'
#' @examples
#' data("rain")
#'
#' ## Fit IDR to data of 185 days using componentwise order on HRES and CTR and
#' ## increasing convex order on perturbed ensemble forecasts (P1, P2, ..., P50)
#'
#' varNames <- c("HRES", "CTR", paste0("P", 1:50))
#' X <- rain[1:185, varNames]
#' y <- rain[1:185, "obs"]
#'
#' ## HRES and CTR are group '1', with componentwise order "comp", perturbed
#' ## forecasts P1, ..., P50 are group '2', with "icx" order
#'
#' groups <- setNames(c(1, 1, rep(2, 50)), varNames)
#' orders <- c("comp" = 1, "icx" = 2)
#'
#' fit <- idr(y = y, X = X, orders = orders, groups = groups)
#' fit
idr <- function(y,
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
                progress = TRUE) {
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
    progress
  )

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
    subsamples     = NULL,
    subsample_size = NULL,
    replace        = NULL,
    settings       = inputs$pars,
    seed           = NULL,
    n_jobs         = 1L,
    show_progress  = inputs$progress
  )

  structure(
    list(
      y = inputs$y,
      X = inputs$X,
      # takes care of duplicate obs as opposed to simply taking the backing cdfs
      cdf = external_ptr$cdf(inputs$X_formatted),
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

validate <- function(y,
                     X,
                     y_observed,
                     weights,
                     decreasing,
                     groups,
                     orders,
                     stoch,
                     pars,
                     progress,
                     seed = NULL) {
  if (!is.vector(y, mode = "numeric")) {
    stop("'y' must be a numeric vector")
  }
  if (length(y) == 0) {
    stop("'y' must contain at least one value")
  }
  if (anyNA(y)) {
    stop("'y' must not contain NAs")
  }
  if (!all(is.finite(y))) {
    stop("'y' must contain only finite values")
  }
  y <- as.double(y)

  if (!is.data.frame(X)) {
    stop("'X' must be a data.frame")
  }
  if (!all(sapply(X, function(col) {
    is.numeric(col) || is.ordered(col)
  }))) {
    stop("'X' must contain numeric or ordered factor variables")
  }
  if (nrow(X) == 0) {
    stop("'X' must have at least 1 row")
  }
  if (nrow(X) != length(y)) {
    stop("length(y) and nrow(X) must match")
  }
  if (anyNA(X)) {
    stop("'X' must not contain NAs")
  }
  if (anyDuplicated(colnames(X))) {
    stop("'X' must not contain duplicated column names")
  }
  # data.matrix encodes ordered factors by their level index; t(X) on a data
  # frame with factor columns would instead go through a character matrix and
  # destroy the values.
  X_formatted <- as.numeric(t(data.matrix(X)))
  if (!all(is.finite(X_formatted))) {
    stop("'X' must contain only finite values")
  }

  if (!is.null(y_observed)) {
    if (length(y_observed) != length(y)) {
      stop("length(y_observed) and length(y) must match")
    }
    if (!(is.logical(y_observed) || is.numeric(y_observed)) ||
          anyNA(y_observed) || !all(y_observed %in% c(0, 1))) {
      stop("'y_observed' must contain only TRUE/1 or FALSE/0 values")
    }
    y_observed <- as.logical(y_observed)
    if (!any(y_observed)) {
      stop("at least one observation must be uncensored")
    }
  }

  if (!is.null(weights)) {
    if (!is.vector(weights, "numeric") || length(weights) != length(y)) {
      stop("'weights' must be a numeric vector as long as 'y'")
    }
    if (anyNA(weights) || !all(is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must contain only finite non-negative values")
    }
    if (!any(weights > 0)) {
      stop("at least one weight must be positive")
    }
    weights <- as.double(weights)
  }

  if (!(identical(decreasing, TRUE) || identical(decreasing, FALSE))) {
    stop("decreasing must be a boolean")
  }

  if (length(orders) > 0 &&
        (is.null(names(orders)) || !all(nzchar(names(orders))))) {
    stop("'orders' must be a named vector")
  }
  if (!all(names(orders) %in% c("comp", "sd", "icx"))) {
    stop("orders must be in 'comp', 'sd', 'icx'")
  }
  if (length(orders) != length(unique(orders))) {
    stop("multiple orders specified for some group(s)")
  }
  if (anyNA(groups)) {
    stop("'groups' must not contain NAs")
  }
  M <- match(colnames(X), names(groups), nomatch = 0)
  if (any(M == 0)) {
    stop("the same variable names must be used in 'groups' and in 'X'")
  }
  if (!all(names(groups) %in% colnames(X))) {
    stop("'groups' contains variable names that are not in 'X'")
  }
  if (!setequal(unname(orders), unname(groups))) {
    stop("different group labels in 'groups' and 'orders'")
  }

  if (!identical(stoch, "sd") && !identical(stoch, "hazard")) {
    stop("only 'sd' or 'hazard' allowed as stochastic order constraints")
  }

  if (isTRUE(progress == 1)) {
    progress <- TRUE
  }
  if (isTRUE(progress == 0)) {
    progress <- FALSE
  }
  if (!isTRUE(progress) && !isFALSE(progress)) {
    stop("'progress' must be TRUE/FALSE or 1/0")
  }

  # Missing options fall back to the documented defaults; this also guarantees
  # the Rust side always receives a complete settings list.
  default_pars <- list(
    verbose = FALSE,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    max_iter = 10000L
  )
  if (is.null(pars)) {
    pars <- default_pars
  } else {
    if (!is.list(pars)) {
      stop("'pars' should be a list of options if provided")
    }
    if (length(pars) > 0 &&
          (is.null(names(pars)) || !all(nzchar(names(pars))))) {
      stop("'pars' options must be named")
    }
    known <- names(pars) %in% names(default_pars)
    if (!all(known)) {
      stop(paste(
        "'pars' option(s)",
        paste(names(pars)[!known], collapse = ", "),
        "unknown"
      ))
    }
    pars <- modifyList(default_pars, pars)
  }
  if (!isTRUE(pars$verbose) && !isFALSE(pars$verbose)) {
    stop("'pars' option 'verbose' should be TRUE or FALSE")
  }
  if (!is.numeric(pars$eps_abs) || length(pars$eps_abs) != 1 ||
        is.na(pars$eps_abs) || pars$eps_abs <= 0.0) {
    stop("'pars' option 'eps_abs' should be a positive number")
  }
  if (!is.numeric(pars$eps_rel) || length(pars$eps_rel) != 1 ||
        is.na(pars$eps_rel) || pars$eps_rel <= 0.0) {
    stop("'pars' option 'eps_rel' should be a positive number")
  }
  if (!is.numeric(pars$max_iter) || length(pars$max_iter) != 1 ||
        is.na(pars$max_iter) || pars$max_iter < 1 ||
        pars$max_iter %% 1 != 0) {
    stop("'pars' option 'max_iter' should be a positive integer")
  }
  pars$max_iter <- as.integer(pars$max_iter)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
      stop("'seed' must be a single non-negative number")
    }
    seed <- as.double(seed)
  }

  list(
    y = y,
    X = X,
    X_formatted = X_formatted,
    y_observed = y_observed,
    weights = weights,
    decreasing = decreasing,
    groups = groups,
    orders = orders,
    stoch = stoch,
    progress = progress,
    pars = pars,
    seed = seed
  )
}

#' Predict method for IDR fits
#'
#' @description Prediction based on IDR model fit.
#'
#' @method predict idrfit
#'
#' @param object IDR fit (object of class \code{"idrfit"}).
#' @param data optional \code{data.frame} containing variables with which to
#'   predict. In-sample predictions are returned if this is omitted. Ordered
#'   factor variables are converted to numeric for computation, so ensure that
#'   the factor levels are identical in \code{data} and the training data for
#'   \code{fit}.
#' @param digits number of decimal places for the predictive CDF. Accepted for
#'   backwards compatibility but currently ignored (a warning is issued once
#'   per session); predictions are returned at full precision.
#' @param interpolation interpolation method for univariate data, ignored at
#'   this time. Only linear is supported for single variate, multivariate uses
#'   midpoint.
#' @param ... included for generic function consistency.
#'
#' @details If the variables \code{x = data[j,]} for which predictions are
#' desired are already contained in the training dataset \code{X} for the fit,
#' \code{predict.idrfit} returns the corresponding in-sample prediction.
#' Otherwise monotonicity is used to derive upper and lower bounds for the
#' predictive CDF, and the predictive CDF is a pointwise average of these
#' bounds. For univariate IDR with a numeric covariate, the predictive CDF is
#' computed by linear interpolation. Otherwise, or if
#' \code{interpolation != "linear"}, midpoint interpolation is used, i.e.
#' default weights of \code{0.5} for both the lower and the upper bound.
#'
#' If the lower and the upper bound on the predictive cdf are far apart (or
#' trivial, i.e. constant 0 or constant 1), this indicates that the prediction
#' based on \code{x} is uncertain because either the training dataset is too
#' small or only few similar variable combinations as in \code{x} have been
#' observed in the training data. However, \emph{the bounds on the predictive
#' CDF are not prediction intervals and should not be interpreted as such. They
#' only indicate the uncertainty of out-of-sample predictions for which the
#' variables are not contained in the training data.}
#'
#' If the new variables \code{x} are greater than all \code{X[i, ]} in the
#' selected order(s), the lower bound on the cdf is trivial (constant 0) and the
#' upper bound is taken as predictive cdf. The upper bound on the cdf is trivial
#' (constant 1) if \code{x} is smaller than all \code{X[i, ]}. If \code{x} is
#' not comparable to any row of \code{X} in the given order, a prediction based
#' on the training data is not possible. In that case, the default forecast is
#' the empirical distribution of \code{y} in the training data.
#'
#' @return A list with the cdf jump points and the values at those jump points
#'   for each covariate.
#'
#' \item{\code{points}}{the points where the predictive CDF has jumps.}
#'
#' \item{\code{cdf}}{the estimated CDF evaluated at the \code{points}.}
#'
#' @export
#' @importFrom stats predict
#' @importFrom stats approx
#'
#' @seealso
#' \code{\link{idr}} to fit IDR to training data.
#'
#' \code{\link{cdf}}, \code{\link{qpred}} to evaluate the CDF or quantile
#' function of IDR predictions.
#'
#' \code{\link{bscore}}, \code{\link{qscore}}, \code{\link{crps}},
#' \code{\link{pit}} to compute Brier scores, quantile scores, the CRPS and the
#' PIT of IDR predictions.
#'
#' \code{\link[isodistrreg:plot.idr]{plot}} to plot IDR predictive CDFs.
#'
#' @examples
#' data("rain")
#'
#' ## Fit IDR to data of 185 days using componentwise order on HRES and CTR and
#' ## increasing convex order on perturbed ensemble forecasts (P1, P2, ..., P50)
#'
#' varNames <- c("HRES", "CTR", paste0("P", 1:50))
#' X <- rain[1:185, varNames]
#' y <- rain[1:185, "obs"]
#'
#' ## HRES and CTR are group '1', with componentwise order "comp", perturbed
#' ## forecasts P1, ..., P50 are group '2', with "icx" order
#'
#' groups <- setNames(c(1, 1, rep(2, 50)), varNames)
#' orders <- c("comp" = 1, "icx" = 2)
#'
#' fit <- idr(y = y, X = X, orders = orders, groups = groups)
#'
#' ## Predict for day 186
#' predict(fit, data = rain[186, varNames])
predict.idrfit <- function(object,
                           data = NULL,
                           digits = NULL,
                           interpolation = NULL,
                           ...) {
  if (is.null(data)) {
    # In-sample predictions, as documented.
    data <- object$X
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (!all(sapply(data, function(col) {
    is.numeric(col) || is.ordered(col)
  }))) {
    stop("'data' must contain numeric or ordered factor variables")
  }
  if (nrow(data) == 0) {
    stop("'data' must have at least 1 row")
  }
  covariate_dimension <- object$external_ptr$dimension()
  if (ncol(data) != covariate_dimension) {
    if (covariate_dimension == 1) {
      message <- "ncol(data) must be 1 for a total order fit"
    } else {
      message <- "ncol(data) must match the covariate dimension"
    }
    stop(message)
  }
  fit_names <- colnames(object$X)
  if (!is.null(fit_names) && !identical(colnames(data), fit_names)) {
    if (setequal(colnames(data), fit_names)) {
      # Same variables in a different order: match by name, not position.
      data <- data[, fit_names, drop = FALSE]
    } else {
      stop("'data' must contain the same variables as the training data")
    }
  }
  if (anyNA(data)) {
    stop("'data' must not contain NAs")
  }
  # See validate(): data.matrix keeps ordered factors usable as level codes.
  new_covariates <- as.numeric(t(data.matrix(data)))
  if (anyNA(new_covariates)) {
    stop("'data' must not contain NAs")
  }

  if (!is.null(digits)) {
    # TODO: Implement?
    warn_once(
      "digits",
      "'digits' parameter is ignored, last warning this session"
    )
  }
  if (!is.null(interpolation)) {
    # TODO: Implement?
    warn_once(
      "interpolation",
      "'interpolation' parameter is ignored, last warning this session"
    )
  }

  preds <- list(
    points = object$response_unique,
    cdf = object$external_ptr$cdf(new_covariates),
    # these two (undocumented) values are used for e.g., quantile predictions
    predict_covariates = new_covariates,
    external_ptr = object$external_ptr
  )

  structure(preds, class = "idr")
}

#' @export
print.idr <- function(x, ...) {
  nr <- NROW(x$cdf)
  nc <- NCOL(x$cdf)

  cat("IDR predictions: list(\n")
  cat(sprintf("  points = vector of %d jump points,\n", nc))
  cdf_line <- "  cdf = matrix of %d prediction(s) at %d jump point(s)\n"
  cat(sprintf(cdf_line, nr, nc))
  cat(")\n")

  cat("\n")

  if (nc > 0L) {
    k <- min(10L, nc)
    vals <- x$points[1:k]
    # Round to two decimals and display with exactly two digits
    disp <- formatC(round(vals, 2), format = "f", digits = 2)
    cat("CDF jump points: ")
    cat(paste(disp, collapse = " "))
    if (nc > k) cat(" ...")
    cat("\n")
  }

  if (nr > 0L && nc > 0L) {
    k <- min(10L, nc)
    vals <- x$cdf[1L, seq_len(k), drop = TRUE]
    # Round to two decimals and display with exactly two digits
    disp <- formatC(round(vals, 2), format = "f", digits = 2)
    cat("CDF values:\n")
    cat("[1] ", paste(disp, collapse = " "), sep = "")
    if (nc > k) cat(" ...")
    cat("\n")
  } else {
    cat("(empty)\n")
  }
  if (nr > 1L) {
    cat("[2] ...\n")
  }

  invisible(x)
}


#' @export
print.idrfit <- function(x, ...) {
  dimension <- x$external_ptr$dimension()
  kind <- if (dimension == 1) {
    "total order"
  } else {
    "partial order"
  }

  cat(paste("IDR fit with", kind, "\n"))

  nr_thresholds <- length(x$response_unique)
  cat(paste("CDFs estimated:", length(x$cdf) / nr_thresholds, "\n"))
  cat(paste("Thresholds for estimation:", nr_thresholds, "\n"))

  if (kind == "partial order") {
    # Print diagnostics info
    prec <- signif(x$diagnostic$precision, 2)
    cat(paste("CDF estimation error:", prec, "\n"))
    conv <- signif(x$diagnostic$convergence_fraction, 4) * 100
    cat(paste0(
      "Converged before hitting max iterations: ",
      conv,
      "% of thresholds\n"
    ))
  }

  invisible(x)
}

#' Isotonic mean regression
#'
#' @description Computes isotonic mean regression for numeric responses. When
#'   covariates are supplied they determine the ordering; when omitted the
#'   responses are assumed pre-sorted (regression on the index). When weights
#'   are omitted every observation receives weight 1.
#'
#' @param y numeric vector of response values.
#' @param X numeric vector of covariate values, or \code{NULL} if responses are
#'   pre-sorted.
#' @param weights numeric vector of finite, non-negative weights, at least one
#'   of which must be positive, or \code{NULL} for equal weights.
#' @param decreasing whether the fit is decreasing in the covariate (default
#'   \code{FALSE} is increasing, \code{TRUE} is decreasing).
#'
#' @return Numeric vector of isotonic fitted means, one per observation.
#'
#' @examples
#' isotonic_regression(c(2, 3, 1, 4, 5), X = 1:5)
#' isotonic_regression(c(3, 2, 4, 1), X = 1:4, weights = c(1, 2, 1, 1))
#' isotonic_regression(sort(c(3, 1, 2, 5)))
#' isotonic_regression(sort(c(2, 1, 3)), weights = c(1, 2, 1))
#'
#' @export
isotonic_regression <- function(y, X = NULL, weights = NULL,
                                decreasing = FALSE) {
  if (!is.vector(y, mode = "numeric")) {
    stop("'y' must be a numeric vector")
  }
  if (anyNA(y) || !all(is.finite(y))) {
    stop("'y' must contain only finite values")
  }
  y <- as.double(y)

  if (!is.null(X)) {
    if (!is.vector(X, mode = "numeric")) {
      stop("'X' must be a numeric vector")
    }
    if (length(X) != length(y)) {
      stop("'X' and 'y' must have equal length")
    }
    if (anyNA(X) || !all(is.finite(X))) {
      stop("'X' must contain only finite values")
    }
    X <- as.double(X)
  }

  if (!is.null(weights)) {
    if (!is.vector(weights, mode = "numeric")) {
      stop("'weights' must be a numeric vector")
    }
    if (length(weights) != length(y)) {
      stop("'weights' and 'y' must have equal length")
    }
    if (anyNA(weights) || !all(is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must contain only finite non-negative values")
    }
    if (length(weights) > 0 && !any(weights > 0)) {
      stop("at least one weight must be positive")
    }
    weights <- as.double(weights)
  }

  if (!isTRUE(decreasing) && !isFALSE(decreasing)) {
    stop("'decreasing' must be TRUE or FALSE")
  }

  isotonic_regression_impl(y, X, weights, decreasing)
}

#' Helper to warn once per session
warn_once <- local({
  seen <- new.env(parent = emptyenv())
  #' Warn only once per session
  #'
  #' @description Filters warnings, letting through only the first of each.
  #'
  #' @param key How to identify the warning.
  #' @param msg Message text.
  #' @param call. Logical; if \code{TRUE}, the call is included in the warning
  #'   message.
  #' @param immediate. Logical; if \code{TRUE}, the warning is issued
  #'   immediately rather than being deferred.
  function(key, msg, call. = FALSE, immediate. = FALSE) {
    if (!exists(key, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      warning(msg, call. = call., immediate. = immediate.)
    }
    invisible(NULL)
  }
})
