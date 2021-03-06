% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bagging.R
\name{idrbag}
\alias{idrbag}
\title{Compute IDR predictions with (su)bagging}
\usage{
idrbag(y, X, groups = setNames(rep(1, ncol(X)), colnames(X)), orders =
  c("comp" = 1), stoch = "sd", pars = osqpSettings(verbose = FALSE, eps_abs =
  1e-5, eps_rel = 1e-5, max_iter = 10000L), progress = TRUE, newdata, 
  digits = 3, interpolation = "linear", b, p, replace = FALSE, grid = NULL)
}
\arguments{
\item{y}{numeric vector (the response variable).}

\item{X}{data frame of numeric or ordered factor variables (the regression
covariates).}

\item{groups}{named vector of length \code{ncol(X)} denoting groups of
variables that are to be ordered with the same order (see 'Details'). Only
relevant if \code{X} contains more than one variable. The same names as in
\code{X} should be used.}

\item{orders}{named vector giving for each group in \code{groups} the order
that will be applied to this group. Only relevant if \code{X} contains more
than one variable. The names of \code{orders} give the order, the entries
give the group labels. Available options: \code{"comp"} for componentwise
order, \code{"sd"} for stochastic dominance, \code{"icx"} for increasing
convex order (see 'Details). Default is \code{"comp"} for all variables.
The \code{"sd"} and \code{"icx"} orders can only be used with numeric
variables, but not with ordered factors.}

\item{stoch}{stochastic order constraint used for estimation. Default is
\code{"sd"} for first order stochastic dominance. Use \code{"hazard"} for
hazard rate order (experimental).}

\item{pars}{parameters for quadratic programming optimization (only relevant
if \code{X} has more than one column), set using
\code{\link[osqp]{osqpSettings}}.}

\item{progress}{display progressbar (\code{TRUE}, \code{FALSE} or \code{1},
\code{0})?}

\item{newdata}{\code{data.frame} containing variables with which to
predict. Ordered factor variables are converted to numeric for computation,
so ensure that the factor levels are identical in \code{newdata} and in
\code{X}.}

\item{digits}{number of decimal places for the predictive CDF.}

\item{interpolation}{interpolation method for univariate data. Default is 
\code{"linear"}. Any other argument will select midpoint interpolation (see 
'Details' in \code{\link{predict.idrfit}}). Has no effect for multivariate
IDR.}

\item{b}{number of (su)bagging samples.}

\item{p}{size of (su)bagging samples relative to training data.}

\item{replace}{draw samples with (\code{TRUE}, \code{1}) or without
(\code{FALSE}, \code{0}) replacement?}

\item{grid}{grid on which the predictive CDFs are evaluated. Default are
the unique values of \code{y}.}
}
\value{
A list of predictions, see \code{\link{predict.idrfit}}.
}
\description{
Computes IDR predictions with bootstrap aggregating (bagging)
or subsample aggregation (subagging).
}
\details{
This function draws \code{b} times a random subsample of size
\code{ceiling(nrow(X)*p)}) from the training data, fits IDR to each
subsample, computes predictions for the new data supplied in \code{newdata},
and averages the predictions derived from the \code{b} subsamples. There are
no default values for \code{b} and \code{p}.
}
