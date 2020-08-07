#' Prepare data for IDR modeling with given orders
#' 
#' @usage
#' prepareData(X, groups, orders)
#' 
#' @param X \code{data.frame} of covariates.
#' @param groups vector of groups.
#' @param orders named vector of orders for groups.
#' 
#' @keywords internal
prepareData <- function(X, groups, orders) {
  grNames <- names(groups)
  ordNames <- names(orders)
  for (group in unique(groups)) {
    M <- grNames[groups == group]
    if (length(M) > 1) {
      if (ordNames[orders == group] == "comp") 
        next
      if (!all(sapply(X[, M], is.numeric))) 
        stop("grouping is only possible for numeric variables")
      tmp <- apply(X[, M], 1, sort, decreasing = TRUE)
      if (ordNames[orders == group] == "sd") 
        X[, M] <- t(tmp)
      else
        X[, M] <- t(apply(tmp, 2, cumsum))
    }
  }
  X
}

#' Fit IDR to training data
#'
#' @description Fits isotonic distributional regression (IDR) to a training
#'   dataset.
#'
#' @usage idr(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)), orders =
#'   c("comp" = 1), stoch = "sd", pars = osqpSettings(verbose = FALSE, eps_abs =
#'   1e-5, eps_rel = 1e-5, max_iter = 10000L), progress = TRUE)
#'
#' @param y numeric vector (the response variable).
#' @param X data frame of numeric or ordered factor variables (the regression
#'   covariates).
#' @param groups named vector of length \code{ncol(X)} denoting groups of
#'   variables that are to be ordered with the same order (see 'Details'). Only
#'   relevant if \code{X} contains more than one variable. The same names as in
#'   \code{X} should be used.
#' @param orders named vector giving for each group in \code{groups} the order
#'   that will be applied to this group. Only relevant if \code{X} contains more
#'   than one variable. The names of \code{orders} give the order, the entries
#'   give the group labels Available options: \code{"comp"} for componentwise
#'   order, \code{"sd"} for stochastic dominance, \code{"icx"} for increasing
#'   convex order (see 'Details). Default is \code{"comp"} for all variables.
#'   The \code{"sd"} and \code{"icx"} orders can only be used with numeric
#'   variables, but not with ordered factors.
#' @param stoch stochastic order constraint used for estimation. Default is
#'   \code{"sd"} for first order stochastic dominance. Use \code{"hazard"} for
#'   hazard rate order (under development, only available for one-dimensional
#'   \code{X}).
#' @param pars parameters for quadratic programming optimization (only relevant
#'   if \code{X} has more than one column), set using
#'   \code{\link[osqp]{osqpSettings}}.
#' @param progress display progressbar?
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
#'   than one variable, the CDFs are estimated by calling
#'   \code{\link[osqp]{solve_osqp}} from the package \pkg{osqp}
#'   \code{length(unique(y))} times. This might take a while if the training
#'   dataset is large.
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
#'   \item{\code{X}}{data frame of all distinct covariate combinations used for
#'   the fit.}
#'
#'   \item{\code{y}}{list of all observed responses in the training data for
#'   given covariate combinations in \code{X}.}
#'
#'   \item{\code{cdf}}{matrix containing the estimated CDFs, one CDF per row,
#'   evaluated at \code{thresholds} (see next point). The CDF in the \code{i}th
#'   row corredponds to the estimated conditional distribution of the response
#'   given the covariates values in \code{X[i,]}.}
#'
#'   \item{\code{thresholds}}{the thresholds at which the CDFs in \code{cdf} are
#'   evaluated. The entries in \code{cdf[,j]} are the conditional CDFs evaluated
#'   at \code{thresholds[j]}.}
#'
#'   \item{\code{groups}, \code{orders}}{ the groups and orders used for
#'   estimation.}
#'
#'   \item{\code{diagnostic}}{list giving a bound on the precision of the CDF
#'   estimation (the maximal downwards-step in the CDF that has been detected)
#'   and the fraction of CDF estimations that were stopped at the iteration
#'   limit \code{max_iter}. Decrease the parameters \code{eps_abs} and/or
#'   \code{eps_rel} or increase \code{max_iter} in \code{pars} to improve the
#'   precision. See \code{\link[osqp]{osqpSettings}} for more optimization
#'   parameters.}
#'
#'   \item{\code{indices}}{ the row indices of the covariates in \code{X} in the
#'   original training dataset (used for in-sample predictions with
#'   \code{\link{predict.idrfit}}).}
#'
#'   \item{\code{constraints}}{ (in multivariate IDR, \code{NULL} otherwise)
#'   matrices giving the order constraints for optimization. Used in
#'   \code{\link{predict.idrfit}}.}
#'
#'
#' @note The function \code{idr} is only intended for fitting IDR model for a
#'   training dataset and storing the results for further processing, but not
#'   for prediction or evaluation, which is done using the output of
#'   \code{\link{predict.idrfit}}.
#'
#' @seealso The S3 method \code{\link{predict.idrfit}} for predictions based on
#'   an IDR fit.
#'
#' @export
#' @importFrom osqp osqp
#' @importFrom stats setNames
#'
#' @references Henzi, A., Moesching, A., & Duembgen, L. (2020). Accelerating the
#' pool-adjacent-violators algorithm for isotonic distributional regression.
#' arXiv preprint arXiv:2006.05527.
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020).
#' OSQP: An operator splitting solver for quadratic programs. Mathematical
#' Programming Computation, 1-36.
#'
#' Bartolomeo Stellato, Goran Banjac, Paul Goulart and Stephen Boyd (2019).
#' osqp: Quadratic Programming Solver using the 'OSQP' Library. R package
#' version 0.6.0.3. https://CRAN.R-project.org/package=osqp
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
idr <- function(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)),
  orders = c("comp" = 1), stoch = "sd", pars = osqpSettings(verbose = FALSE,
  eps_abs = 1e-5, eps_rel = 1e-5, max_iter = 10000L), progress = TRUE) {
    
  # Check input
  if (!is.vector(y, mode = "numeric")) 
    stop("'y' must be a numeric vector")
  if (!is.data.frame(X)) 
    stop("'X' must be a data.frame")
  if (!all(sapply(X, function(col) is.numeric(col) || is.ordered(col)))) 
    stop("'X' must contain numeric or ordered factor variables")
  if (nrow(X) <= 1) 
    stop("'X' must have more than 1 rows")
  if (nrow(X) != length(y)) 
    stop("length(y) and nrow(X) must match")
  if (anyNA(X) || anyNA(y)) 
    stop("'X' and 'y' must not contain NAs")
  if (!all(names(orders) %in% c("comp", "sd", "icx"))) 
    stop("orders must be in 'comp', 'sd', 'icx'")
  if (length(orders) != length(unique(orders)))
    stop("multiple orders specified for some group(s)")
  M <- match(colnames(X), names(groups), nomatch = 0)
  if (any(M == 0))
    stop("the same variable names must be used in 'groups' and in 'X'")
  if (!identical(sort(unique(unname(orders))), sort(unique(groups)))) 
    stop("different group labels in 'groups' and 'orders")
  thresholds <- sort(unique(y))
  if ((nThr <- length(thresholds)) == 1) 
    stop("'y' must contain more than 1 distinct value")
  if (ncol(X) > 1 & !identical(stoch, "sd"))
    stop("only first order stochastic dominance for multivariate X available")
  if (!identical(stoch, "sd") & ! identical(stoch, "hazard"))
    stop("only 'sd' or 'hazard' allowed as stochastic order constraints")
  if (!isTRUE(progress) & !isFALSE(progress))
    stop("'progress' must be TRUE or FALSE")
  X <- prepareData(X, groups, orders)
  
  # Aggregate input data
  nVar <- ncol(X)
  oldNames <- names(X)
  names(X) <- NULL
  x <- data.frame(y = y, ind = seq_along(y))
  X <- stats::aggregate(x = x, by = X, FUN = identity, simplify = FALSE)
  cpY <- X[["y"]]
  indices <- X[["ind"]]
  X <- X[, 1:nVar, drop = FALSE]
  names(X) <- oldNames
  weights <- sapply(indices, length)
  
  if (nVar == 1) {
    # One-dimensional IDR using PAVA
    constr <- NULL
    diagnostic <- list(precision = 0, convergence = 0)
    if (stoch == "sd") {
      cdf <- isoCdf_sequential(
        w = weights,
        W = rep(1, length(y)),
        Y = sort(y),
        posY = rep.int(seq_along(indices),lengths(indices))[order(unlist(cpY))],
        y = thresholds
      )$CDF
    } else {
      cdf <- idrHazardCpp(
        w = weights,
        W = rep(1, length(y)),
        Y = sort(y),
        posY = rep.int(seq_along(indices),lengths(indices))[order(unlist(cpY))],
        y = thresholds
      )
    }
  } else {
    # Multivariate IDR using osqp
    constr <- compOrd(X)
    N <- nrow(X)
    cdf <- matrix(ncol = nThr - 1, nrow = N)
    A <- trReduc(constr$paths, N)
    nConstr <- nrow(A)
    l <- rep(0, nConstr)
    A <- Matrix::sparseMatrix(i = rep(seq_len(nConstr), 2), j = as.vector(A), 
      x = rep(c(1, -1), each = nConstr), dims = c(nConstr, N))
    P <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = weights)
    i <- 1
    I <- nThr - 1
    conv <- vector("logical", I)
    
    q <- - weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
    qp <- osqp::osqp(P = P, q = q, A = A, l = l, pars = pars)
    sol <- qp$Solve()
    cdf[, 1] <- pmin(1, pmax(0, sol$x))
    conv[1] <- identical(sol$info$status, "maximum iterations reached")
    
    if (I > 1) {
      if (progress) {
        cat("Estimating cdf...\n")
        pb <- utils::txtProgressBar(style = 1)
        for (i in 2:I) {
          utils::setTxtProgressBar(pb, i/I)
          qp$WarmStart(x = cdf[, i - 1L])
          q <-  -weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
          qp$Update(q = q)
          sol <- qp$Solve()
          cdf[, i] <- pmin(1, pmax(0, sol$x))
          conv[i] <- identical(sol$info$status, "maximum iterations reached")
        }
        close(pb)
        cat("\n")
      } else {
        for (i in 2:I) {
          qp$WarmStart(x = cdf[, i - 1L])
          q <-  -weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
          qp$Update(q = q)
          sol <- qp$Solve()
          cdf[, i] <- pmin(1, pmax(0, sol$x))
          conv[i] <- identical(sol$info$status, "maximum iterations reached")
        }
      }
      diagnostic <- list(
        precision = ifelse(I > 1, abs(min(diff(t(cdf)))), 0),
        convergence = mean(conv)
      )
    }
  }
  
  # Apply pava to estimated CDF to ensure monotonicity
  if (nVar > 1) cdf <- cbind(pavaCorrect(cdf), 1)
  
  structure(list(X = X, y = cpY, cdf = cdf, thresholds = thresholds, 
    groups = groups, orders = orders, diagnostic = diagnostic,
    indices = indices, constraints = constr), class = "idrfit")
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
#' @param digits number of decimal places for the predictive CDF.
#' @param interpolation interpolation method for univariate data. Default is 
#'   \code{"linear"}. Any other argument will select midpoint interpolation (see 
#'   'Details'). Has no effect for multivariate IDR.
#' @param asplitAvail use \code{\link[base]{asplit}} for splitting arrays
#'   (default is \code{TRUE}). Set to \code{FALSE} for R Versions < 3.6, where
#'   \code{asplit} is not available.
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
#' @return A list of predictions. Each prediction is a \code{data.frame}
#'   containing the following variables:
#'
#' \item{\code{points}}{the points where the predictive CDF has jumps.}
#'
#' \item{\code{cdf}}{the estimated CDF evaluated at the \code{points}.}
#'
#' \item{\code{lower}, \code{upper}}{ (only for out-of-sample predictions)
#'   bounds for the estimated CDF, see 'Details' above.}
#'
#' The output has the attribute \code{incomparables}, which gives the indices
#' of all predictions for which the climatological forecast is returned because
#' the forecast variables are not comparable to the training data.
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
predict.idrfit <- function(object, data = NULL, digits = 3,
    interpolation = "linear", asplitAvail = TRUE, ...) {
  cdf <- object$cdf
  thresholds <- object$thresholds
  
  if (is.null(data)) {
    # In-sample predictions
    indices <- object$indices
    preds <- structure(vector("list", length(unlist(indices))),
        class = c("idr"))
    for (i in seq_along(indices)) {
      edf <- round(cdf[i, ], digits)
      sel <- c(edf[1] > 0, diff(edf) > 0)
      tmp <- data.frame(points = thresholds[sel], cdf = edf[sel])
      for (j in indices[[i]]) {
        preds[[j]] <- tmp
      }
    }
    return(preds)
  }
  
  # Out-of-sample predictions
  if (!is.data.frame(data)) 
    stop("'data' must be a data.frame")
  X <- object$X
  M <- match(colnames(X), colnames(data), nomatch = 0)
  if (any(M == 0)) 
    stop("some variables of the idr fit are missing in 'data'")
  data <- prepareData(data[, M, drop = FALSE], groups = object$groups, 
    orders = object$orders)
  nVar <- ncol(data)
  nx <- nrow(data)
  
  # Prediction method for univariate IDR
  if (nVar == 1) {
    X <- X[[1]]
    x <- data[[1]]
    fct <- is.ordered(X)
    if (fct) {
      X <- as.integer(X)
      x <- as.integer(x)
    }
    smaller <- findInterval(x, X)
    smaller[smaller == 0] <- 1
    wg <- stats::approx(x = X, y = seq_along(X), xout = x, yleft = 1, 
      yright = length(X), rule = 2)$y - seq_along(X)[smaller]
    greater <- smaller + as.integer(wg > 0)
    m <- length(thresholds)
    if (identical(interpolation, "linear") & !fct) {
      ws <- 1 - wg
    } else {
      wg <- ws <- rep(0.5, length(x))     
    }
    preds <- Map(
      function(l, u, ws, wg) {
        ind <- (c(0, l[-m]) < l) | (c(0, u[-m]) < u)
        l <- l[ind]
        u <- u[ind]
        cdf <- round(l * wg + u * ws, digits)
        data.frame(
          points = thresholds[ind],
          lower = l,
          cdf = cdf,
          upper = u
        )
      },
      l = splitArr(round(cdf[greater, , drop = FALSE], digits), 1, asplitAvail),
      u = splitArr(round(cdf[smaller, , drop = FALSE], digits), 1, asplitAvail),
      ws = ws,
      wg = wg
    )
    return(structure(preds, class = "idr", incomparables = integer(0)))
  }
  
  # Prediction Method for multivariate IDR
  preds <- structure(vector("list", nx), class = c("idr"))
  nPoints <- neighborPoints(x = data.matrix(data), X = data.matrix(X), 
    orderX = object$constraints, asplitAvail = asplitAvail)
  smaller <- nPoints$smaller
  greater <- nPoints$greater
  
  # Climatological forecast for incomparable variables
  incomparables <- sapply(smaller, length) + sapply(greater, length) == 0
  if (any(incomparables)) {
    y <- unlist(object$y)
    edf <- round(stats::ecdf(y)(thresholds), digits)
    sel <- edf > 0
    edf <- edf[sel]
    points <- thresholds[sel]
    upr <- which.max(edf == 1)
    if (upr < length(edf)) {
      points <- points[-((upr + 1):length(edf))]
      edf <- edf[-((upr + 1):length(edf))]
    }
    dat <- data.frame(points = points, lower = edf, cdf = edf, upper = edf)
    for (i in which(incomparables)) preds[[i]] <- dat
  }
  
  # Predictions for comparable variables
  for (i in which(!incomparables)) {
    # Case distinction: Existence of lower and/or upper bound on CDF
    if (length(smaller[[i]]) > 0 && length(greater[[i]]) == 0) {
      upper <- round(apply(cdf[smaller[[i]], , drop = FALSE], 2, min), digits)
      sel <- c(upper[1] != 0, diff(upper) != 0)
      upper <- estimCdf <- upper[sel]
      lower <- rep(0, length(upper))
    } else if (length(smaller[[i]]) == 0 && length(greater[[i]]) > 0) {
      lower <- round(apply(cdf[greater[[i]], , drop = FALSE], 2, max), digits)
      sel <- c(lower[1] != 0, diff(lower) != 0)
      lower <- estimCdf <- lower[sel]
      upper <- rep(1, length(lower))
    } else {
      lower <- round(apply(cdf[greater[[i]], , drop = FALSE], 2, max), digits)
      upper <- round(apply(cdf[smaller[[i]], , drop = FALSE], 2, min), digits)
      
      sel <- c(lower[1] != 0, diff(lower) != 0) |
        c(upper[1] != 0, diff(upper) != 0)
      lower <- lower[sel]
      upper <- upper[sel]
      estimCdf <- round((lower + upper) / 2, digits)
    }
    preds[[i]] <- data.frame(points = thresholds[sel], lower = lower, 
      cdf = estimCdf, upper = upper)
  }
  attr(preds, "incomparables") <- which(incomparables)
  preds
}

#' @export
print.idr <- function(x, ...) {
  # Print only head of first prediction
  cat(paste("IDR predictions: list of", length(x), "prediction(s)\n\n"))
  print(list(utils::head(x[[1]])))
  cat("...")
  invisible(x)
}

#' @export
print.idrfit <- function(x, ...) {
  # Print diagnostic information
  cat("IDR fit: \n")
  cat(paste("CDFs estimated:", nrow(x$cdf), "\n"))
  cat(paste("Thresholds for estimation:", ncol(x$cdf), "\n"))
  prec <- signif(x$diagnostic$precision, 2)
  conv <- signif(x$diagnostic$convergence, 4) * 100
  cat(paste("CDF estimation precision:", prec, "\n"))
  cat(paste0("Estimations stopped after max_iter iterations: ", conv, "% \n"))
  invisible(x)
}

#' @importFrom osqp osqpSettings
#' @export
osqp::osqpSettings