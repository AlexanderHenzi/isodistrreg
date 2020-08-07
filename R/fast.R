#' Fast version of IDR (experimental)
#'
#' @description Fits isotonic distributional regression (IDR) to a training
#' dataset. Is faster than \code{\link{idr}} for univariate data. \emph{This 
#' function is still under development.}
#'
#' @usage idrF(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)),
#' orders = c("comp" = 1), pars = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
#'   eps_rel = 1e-5, max_iter = 10000L))
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
#'   give the group labels Available options: \code{"comp"} for
#'   componentwise order, \code{"sd"} for stochastic dominance, \code{"icx"} for
#'   increasing convex order (see 'Details). Default is \code{"comp"} for
#'   all variables. The \code{"sd"} and \code{"icx"} orders can only be used
#'   with numeric variables, but not with ordered factors.
#' @param pars parameters for quadratic programming optimization (only relevant
#'   if \code{X} has more than one column), set using
#'   \code{\link[osqp]{osqpSettings}}.
#'
#' @details This function computes the isotonic distributional regression (IDR)
#' of a response \emph{y} on on one or more covariates \emph{X}. IDR estimates
#' the cumulative distribution function (CDF) of \emph{y} conditional on
#' \emph{X} by monotone regression, assuming that \emph{y} is more likely to
#' take higher values, as \emph{X} increases. Formally, IDR assumes that the
#' conditional CDF \eqn{F_{y | X = x}(z)} at each fixed threshold \emph{z}
#' decreases, as \emph{x} increases, or equivalently, that the exceedance
#' probabilities for any threshold \code{z} \eqn{P(y > z | X = x)} increase with
#' \emph{x}.
#'
#' The conditional CDFs are estimated at each threshold in \code{unique(y)}.
#' This is the set where the CDFs may have jumps. If \code{X} contains more than
#' one variable, the CDFs are estimated by calling
#' \code{\link[osqp]{solve_osqp}} from the package \pkg{osqp}
#' \code{length(unique(y))} times. This might take a while if the training
#' dataset is large.
#'
#' Use the argument \code{groups} to group \emph{exchangeable} covariates.
#' Exchangeable covariates are indistinguishable except from the order in which
#' they are labelled (e.g. ensemble weather forecasts, repeated measurements
#' under the same measurement conditions).
#'
#' The following orders are available to perform the monotone regression in IDR:
#' \itemize{
#' \item Componentwise order (\code{"comp"}): A covariate vector \code{x1} is
#'   greater than \code{x2} if \code{x1[i] >= x2[i]} holds for all components
#'   \code{i}. This is the \emph{standard order used in multivariate monotone
#'   regression} and \emph{should not be used for exchangeable variables (e.g.
#'   perturbed ensemble forecasts)}.
#' \item Stochastic dominance (\code{"sd"}): \code{x1} is greater than \code{x2}
#'   in the stochastic order, if the (empirical) distribution of the elements of
#'   \code{x1} is greater than the distribution of the elements of \code{x2} (in
#'   first order) stochastic dominance. The \code{"sd"} order is invariant under
#'   permutations of the grouped variables and therefore \emph{suitable for
#'   exchangeable covariables}.
#' \item Increasing convex order (\code{"icx"}): The \code{"icx"} order can
#'   be used for groups of exchangeable variables. It should be used if the
#'   variables have increasing variability, when their mean increases (e.g.
#'   precipitation forecasts or other variables with right-skewed
#'   distributions). More precisely, \code{"icx"} uses the increasing convex
#'   stochastic order on the empirical distributions of the grouped variables.
#' }
#'
#' @return An object of class \code{"idrfit"} containing the following
#' components:
#'
#' \item{\code{X}}{data frame of all distinct covariate combinations used for
#'   the fit.}
#'
#' \item{\code{y}}{list of all observed responses in the training data for given
#'   covariate combinations in \code{X}.}
#'
#' \item{\code{cdf}}{matrix containing the estimated CDFs, one CDF per row,
#'   evaluated at \code{thresholds} (see next point). The CDF in the \code{i}th
#'   row corredponds to the estimated conditional distribution of the response
#'   given the covariates values in \code{X[i,]}.}
#'
#' \item{\code{thresholds}}{the thresholds at which the CDFs in \code{cdf} are
#'   evaluated. The entries in \code{cdf[,j]} are the conditional CDFs evaluated
#'   at \code{thresholds[j]}.}
#'
#' \item{\code{groups}, \code{orders}}{ the groups and orders used for
#'   estimation.}
#'
#' \item{\code{diagnostic}}{list giving a bound on the precision of the CDF
#'   estimation (the maximal downwards-step in the CDF that has been detected)
#'   and the fraction of CDF estimations that were stopped at the iteration
#'   limit \code{max_iter}. Decrease the parameters \code{eps_abs} and/or
#'   \code{eps_rel} or increase \code{max_iter} in \code{pars} to improve the
#'   precision. See \code{\link[osqp]{osqpSettings}} for more optimization
#'   parameters.}
#'
#' \item{\code{indices}}{ the row indices of the covariates in \code{X} in the
#'   original training dataset (used for in-sample predictions with
#'   \code{\link{predict.idrfit}}).}
#'
#' \item{\code{constraints}}{ (in multivariate IDR, \code{NULL} otherwise)
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
#' @references
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. y Boyd, S.. "OSQP: An
#' Operator Splitting Solver for Quadratic Programs". ArXiv e-prints. 2017
#' 
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. y Boyd, S.. osqp:
#' Quadratic Programming Solver using the 'OSQP' Library. 2018, R package
#' version 0.5.0. \url{https://CRAN.R-project.org/package=osqp}.
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
idrF <- function(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)),
  orders = c("comp" = 1), pars = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
  eps_rel = 1e-5, max_iter = 10000L)) {
    
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
  weights <- lengths(indices)
  
  if (nVar == 1) {
    # One-dimensional IDR using PAVA
    constr <- NULL
    diagnostic <- list(precision = 0, convergence = 0)
    cdf <- isoCdf_sequential(
      w = weights,
      W = rep(1, length(y)),
      Y = sort(y),
      posY = rep.int(seq_along(indices), lengths(indices))[order(unlist(cpY))],
      y = thresholds
    )$CDF
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
    
    cat("Estimating cdf...\n")
    pb <- utils::txtProgressBar(style = 3)
    q <- - weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
    qp <- osqp::osqp(P = P, q = q, A = A, l = l, pars = pars)
    sol <- qp$Solve()
    cdf[, 1] <- pmin(1, pmax(0, sol$x))
    conv[1] <- identical(sol$info$status, "maximum iterations reached")
    
    if (I > 1) {
     for (i in 2:I) {
      utils::setTxtProgressBar(pb, i/I)
      qp$WarmStart(x = cdf[, i - 1L])
      q <-  -weights * sapply(cpY, FUN = function(x) mean(thresholds[i] >= x))
      qp$Update(q = q)
      sol <- qp$Solve()
      cdf[, i] <- pmin(1, pmax(0, sol$x))
      conv[i] <- identical(sol$info$status, "maximum iterations reached")
     }
    }
    close(pb)
    cat("\n")
    diagnostic <- list(
      precision = ifelse(I > 1, abs(min(diff(t(cdf)))), 0),
      convergence = mean(conv)
    )
  }
  
  # Apply pava to estimated CDF to ensure monotonicity
  if (nVar > 1) {
    cdf <- pavaCorrect(cdf)
    cdf <- cbind(cdf, 1)
  }
  
  structure(list(X = X, y = cpY, cdf = cdf, thresholds = thresholds, 
    groups = groups, orders = orders, diagnostic = diagnostic,
    indices = indices, constraints = constr), class = "idrfit")
}