#------------------------------------------------------------------------------#
# These functions are used in "simulation_example.R" for the computation of
# distributional regression models in the simulation study of

# Alexander Henzi, Johanna F. Ziegel, and Tilmann Gneiting. Isotonic
# Distributional Regression. arXiv e-prints, art. arXiv:1909.03725, Sep 2019.
# URL https://arxiv.org/abs/1909.03725.

# File is available on GitHub: https://github.com/AlexanderHenzi/isodistrreg.

# Last revised 2021/03/21

#------------------------------------------------------------------------------#
# fit distributional regression models

probFit <- function(method, data, ...) {
  # check data format
  if (!is.data.frame(data) || !identical("y", colnames(data)[1]))
    stop("Check input data!")
  
  # dot-arguments
  dotArgs <- list(...)
  
  # regression methods
  if (identical(method, "idr")) {
    y <- data[, 1]
    X <- data[, -1, drop = FALSE]
    idr(y = y, X = X)
  } else if (identical(method, "idrbag")) {
    data[, -1, drop = FALSE]
  } else if (identical(method, "np")) {
    ydat <- data$y
    xdat <- data[, -1, drop = FALSE]
    npcdistbw(
      ydat = ydat,
      xdat = xdat,
      oykertype = dotArgs$oykertype,
      nmulti = dotArgs$nmulti,
      bwtype = bwtype
    )
  } else if (identical(method, "tram")) {
    if (is.ordered(data$y)) {
      Polr(
        y ~ .,
        data = data        
      )
    } else {
      Colr(
        y ~ .,
        data = data
      )
    }
  } else if (identical(method, "sqr")) {
    rq(
      y ~ .,
      data = data,
      tau = dotArgs$tau
    )
  } else if (identical(method, "qrf")) {
    quantile_forest(
      X = matrix(data$x, ncol = 1),
      Y = data$y,
      quantiles = dotArgs$quantiles,
      min.node.size = dotArgs$min.node.size,
      sample.fraction = dotArgs$sample.fraction,
      honesty.fraction = dotArgs$honesty.fraction,
      honesty = dotArgs$honesty
    )
  }
}

#------------------------------------------------------------------------------#
# get predictions from distributional regression model
probPredict <- function(fit, newdata, ...) {
  UseMethod("probPredict")
}

# method for idr
probPredict.idrfit <- function(fit, newdata, ...) {
  # dot-arguments
  dotArgs <- list(...)
  digits <- dotArgs$digits

  # get predictions
  predict(
    object = fit,
    data = newdata,
    digits = digits,
    interpolation = "linear"
  )
}

# method for idr bagging (directly from data.frame)
probPredict.data.frame <- function(fit, newdata, ...) {
  # dot-arguments
  dotArgs <- list(...)
  y <- dotArgs$y
  digits <- dotArgs$digits
  b <- dotArgs$b
  
  # get predictions
  interpolation <- "linear"
  X <- fit
  grd <- sort(unique(y))
  N <- nrow(X)
  preds <- matrix(nrow = nrow(newdata), ncol = length(grd), 0)
  for (i in seq_len(b)) {
    sel <- sample(N, ceiling(N / 2))
    fit <- idr(y = y[sel], X = X[sel, , drop = FALSE])
    prd <- predict(fit, newdata, digits = digits, interpolation = interpolation)
    preds <- preds + cdf(prd, grd)
  }
  preds <- preds / b
  
  # output
  cdf <- lapply(
    asplit(round(preds, digits), 1),
    function(x) data.frame(points = grd, cdf = x)
  )
  cdf <- lapply(
    cdf,
    function(x) x[c(x$cdf[1] > 0, x$cdf[-1] > x$cdf[-nrow(x)]), ]
  )
  structure(unname(cdf), class = c("idr", "idrbag"))
}

# method for np
probPredict.condbandwidth <- function(fit, newdata, ...) {
  # dot-arguments
  dotArgs <- list(...)
  grd <- dotArgs$grd
  digits <- dotArgs$digits
  
  # get predictions
  ngrd <- length(grd)
  nnew <- nrow(newdata)
  preds <- npcdist(
    bws = fit,
    exdat = newdata[rep(seq_len(nnew), each = ngrd), , drop = FALSE],
    eydat = rep(grd, nnew)
  )
  
  # output
  probs <- split(round(preds$condist, digits), rep(seq_len(nnew), each = ngrd))
  cdf <- lapply(probs, function(x) data.frame(points = grd, cdf = x))
  cdf <- lapply(
    cdf,
    function(x) x[c(FALSE, x$cdf[-1] > x$cdf[-nrow(x)]), , drop = FALSE]
  )
  structure(unname(cdf), class = c("idr", "nppred"))
}

# method for tram
probPredict.tram <- function(fit, newdata, ...) {
  # dot-arguments
  dotArgs <- list(...)
  K <- dotArgs$K
  digits <- dotArgs$digits
  
  # get predictions
  preds <- predict(
    object = fit,
    newdata = newdata,
    type = "distribution",
    K = K
  )
  
  # output
  points <- as.numeric(rownames(preds))
  probs <- asplit(round(preds, digits), 2)
  cdf <- lapply(probs, function(x) data.frame(points = points, cdf = x))
  cdf <- lapply(
    cdf, 
    function(x) x[c(x$cdf[1] > 0, x$cdf[-1] > x$cdf[-nrow(x)]), , drop = FALSE]
  )
  structure(unname(cdf), class = c("idr", "trampred"))
}

# method for quantreg
probPredict.rqs <- function(fit, newdata, ...) {
  # dot-arguments
  dotArgs <- list(...)

  # get predictions
  preds <- t(apply(predict(fit, newdata = newdata), 1, sort))
  
  # output
  qgrd <- fit$tau
  qgrd[length(qgrd)] <- 1
  cdf <- lapply(
    asplit(preds, 1), 
    function(x) data.frame(points = x, cdf = qgrd)
  )
  structure(unname(cdf), class = c("idr", "rqpred"))
}

# method for quantile random forest
probPredict.quantile_forest <- function(fit, newdata, ...) {
  # dot-arguments
  dotArgs <- list(...)

  # get predictions
  qgrd <- dotArgs$tau
  preds <- predict(fit, matrix(newdata$x, ncol = 1), quantiles = qgrd)
  
  # output
  qgrd[length(qgrd)] <- 1
  cdf <- lapply(
    asplit(preds, 1), 
    function(x) data.frame(points = x, cdf = qgrd)
  )
  structure(unname(cdf), class = c("idr", "rqpred"))  
}