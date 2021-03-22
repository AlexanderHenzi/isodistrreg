#' Distributional index model (DIM)
#'
#' @description Fits distributional index model with user-specified index
#' function to training dataset. See the examples at the bottom to learn
#' how to specify a distributional single index model.
#' 
#' @param formula object of class \code{formula} that describes the index 
#'   model
#' @param indexfit function that fits the index model to training data. Should
#'   accept arguments \code{formula} and \code{data} and admit a \code{predict}
#'   method. Further arguments in \code{...} are passed to indexfit.
#'   See examples.
#' @param data \code{data.frame} containing the covariates of the index model
#'   and the response variable.
#' @param response name of the response variable in \code{data}.
#' @param pars parameters for quadratic programming optimization (only relevant
#'   for multivariate index functions), set using
#'   \code{\link[osqp]{osqpSettings}}.
#' @param progress display progressbar for fitting idr?
#' @param ... further arguments passed to \code{indexfit}.
#'
#' @details
#' This function fits a distributional index model (DIM) to training data. The
#' DIM assumes that the response is more likely to attain higher values when the
#' values of the index function increases. The index function can be
#' estimated by parametric methods like \code{\link[stats]{lm}} or 
#'  \code{\link[stats]{glm}} or also nonparametrically.
#'  
#' The formal mathematical assumption of the DIM is that the conditional CDFs
#' \eqn{F_{y | g(X) = g(x)}(z)} at each fixed threshold \emph{z} decreases, as
#' \emph{g(x)} increases. Here \code{y} denotes the response, \code{x}, \code{X}
#' are the covariates in \code{data} and \code{g} is the index function
#' estimated by \code{indexfit}.
#' 
#' Estimation is performed in two steps: \code{indexfit} is applied to 
#' \code{data} to estimate the function \code{g}. With this estimate,
#' \code{\link{idr}} is applied with the pseudo-covariates \code{g(x)} and
#' response \code{y}.
#' 
#' @return 
#' Object of class \code{dindexm}: A list containing the index model (first
#' component) and the IDR fit on the pseudo-data with the index as covariate
#' (second component).
#' 
#' @seealso 
#' \code{\link{idr}} for more information on IDR,
#' \code{\link{predict.dindexfit}} for (out-of-sample) predictions based on a
#' model with with \code{dindexm}.
#' 
#' @references
#' Henzi, A., Kleger, G. R., & Ziegel, J. F. (2020). Distributional (Single)
#' Index Models. arXiv preprint arXiv:2006.09219.
#'
#' @export
#'
#' @examples
#' n <- 1000
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#' y <- rnorm(n, 1 - X[, 1] + X[, 2]^2 / 3 - (1 - X[, 3]) * (1 + X[, 3]) / 2)
#' data <- cbind(y = y, as.data.frame(X))
#' 
#' ## data for out-of-sample prediction
#' newX <- data.frame(x1 = rnorm(10), x2 = rnorm(10), x3 = rnorm(10))
#' 
#' ## linear regression model for index
#' model <- dindexm(
#'   formula = y ~ poly(x1, degree = 2) + poly(x2, degree = 2) + 
#'     poly(x3, degree = 2),
#'   indexfit = lm,
#'   response = "y",
#'   data = data
#' )
#' pred <- predict(model, data = newX)
#' 
#' ## plot
#' plot(pred, 1, main = "LM based DIM")
#' grd <- pred[[1]]$points
#' trueCdf <- pnorm(
#'   grd,
#'   1 - newX[1, 1] + newX[1, 2]^2 / 3 - (1 - newX[1, 3]) * (1 + newX[1, 3]) / 2
#' )
#' points(grd, trueCdf, type = "l", col = 2)
dindexm <- function(formula, indexfit, data, response,
  pars = osqpSettings(verbose = FALSE, eps_abs = 1e-5, eps_rel = 1e-5,
  max_iter = 10000L), progress = TRUE, ...) {

  indexFit <- indexfit(formula = formula, data = data, ...)
  X <- predict(indexFit)
  if (is.numeric(X) && !is.matrix(X) && !is.vector(X) && !is.array(X))
    stop("predict method for 'indexfit' must return numeric matrix or vector")
  if (is.array(X) && length(attributes(X)$dimnames) == 1) {
    X <- matrix(c(X), ncol = 1)
  }
  X <- as.data.frame(X)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  idrFit <- idr(y = data[, response, drop = TRUE], X = X, pars = pars,
    progress = progress)
  structure(list(indexFit = indexFit, idrFit = idrFit), class = "dindexfit")
}

#' Predict method for distributional index model (DIM)
#'
#' @description Prediction based on distributional index model fit.
#'
#' @method predict dindexfit
#'
#' @param object DIM fit (object of class \code{"dindexfit"}).
#' @param data optional \code{data.frame} containing variables with which to
#'   predict. In-sample predictions are returned if this is omitted.
#' @param digits number of decimal places for the predictive CDF.
#' @param interpolation interpolation method for univariate index Default is 
#'   \code{"linear"}. Any other argument will select midpoint interpolation (see 
#'   'Details' in \code{\link{predict.idrfit}}). Has no effect for multivariate
#'   index function.
#' @param asplitAvail use \code{\link[base]{asplit}} for splitting arrays
#'   (default is \code{TRUE}). Set to \code{FALSE} for R Versions < 3.6, where
#'   \code{asplit} is not available.
#' @param ... further arguments passed to the index prediction function.
#' 
#' @return
#' A list of predictions, as for \code{\link{predict.idrfit}}.
#' 
#' @export
#' 
#' @seealso 
#' Examples in \code{\link{dindexm}}.
predict.dindexfit <- function(object, data = NULL, digits = 3,
  interpolation = "linear", asplitAvail = TRUE, ...) {
  
  indexFit <- object$indexFit
  idrFit <- object$idrFit
  if (is.null(data)) {
    return(predict(idrFit, digits = digits, interpolation = interpolation,
      asplitAvail = asplitAvail))
  }

  indexPred <- as.data.frame(predict(indexFit, data, ...))
  if (is.array(indexPred)) {
    colnx <- length(attributes(indexPred)$dimnames)
    indexPred <- matrix(c(indexPred), ncol = colnx)
  }
  indexPred <- as.data.frame(indexPred)
  colnames(indexPred) <- paste0("x", seq_len(ncol(indexPred)))
  
  predict(object = idrFit, data = indexPred, digits = digits,
    interpolation = interpolation, asplitAvail = asplitAvail)
}
