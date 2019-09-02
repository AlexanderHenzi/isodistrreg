#' Isotonic distributional regression (IDR)
#' 
#' Isotonic distributional Regression (IDR) is a nonparametric forecast
#' calibration method based on monotone regression.
#' 
#' @section How does it work?:
#' 
#' Link to paper.
#' 
#' @section The \pkg{isodistrreg} package:
#' 
#' To make probabilistic forecasts with IDR,
#' \itemize{
#' \item call \code{\link{idr}(y = y, X = X, ...)}, where \code{y} is the
#'   response variable (e.g. weather variable observations) and \code{X} is a
#'   \code{data.frame} of covariates (e.g. ensemble forecasts).
#' \item use \code{\link[=predict.idrfit]{predict}(fit, data)}, where \code{fit}
#'   is the model fit computed with \code{idr} and \code{data} is the data based
#'   on which you want to make predictions.
#' }
#' The following pre-defined functions are available to evaluate IDR
#' predictions:
#' \itemize{
#' \item \code{\link{cdf}} and \code{\link{qpred}} to compute the cumulative
#' distribution function (CDF) and quantile function of IDR predictions.
#' \item \code{\link{bscore}} and \code{\link{qscore}} to calculate Brier scores
#'   for probability forecasts for threshold exceedance (e.g. probability of
#'   precipitation) and quantile scores (e.g. mean absolute error of median
#'   forecast.)
#' \item \code{\link{crps}} to compute the continuous ranked probability score
#' (CRPS).
#' \item \code{\link{pit}} to compute the probability integral transform (PIT).
#' \item \code{\link[=plot.idr]{plot}} to plot IDR predictive CDFs.
#' }
#' Use the dataset \code{\link{rain}} to test IDR.
#' 
#' @name isodistrreg-package
#' @author
#' 
#' Package: Alexander Henzi 
#' @useDynLib isodistrreg, .registration = TRUE
#' @name isodistrreg
NULL

#' Unload dll when package is unloaded
#' @keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("isodistrreg", libpath)
}