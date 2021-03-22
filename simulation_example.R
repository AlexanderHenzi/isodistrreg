#------------------------------------------------------------------------------#
# This code was used to compute the simulation examples in 

# Alexander Henzi, Johanna F. Ziegel, and Tilmann Gneiting. Isotonic
# Distributional Regression. arXiv e-prints, art. arXiv:1909.03725, Sep 2019.
# URL https://arxiv.org/abs/1909.03725.

# File is available on GitHub: https://github.com/AlexanderHenzi/isodistrreg.

# Computation may take a few minutes depending on your CPU. Reduce "n" (size of
# training data) and "nValid" (size of validation data) to increase the speed.

# Last revised 2021/03/21

#------------------------------------------------------------------------------#
# packages (must be installed)
library(tidyverse)
library(scoringRules)
library(isodistrreg)
library(np)
library(tram)
library(quantreg)
library(splines)
library(grf)

#------------------------------------------------------------------------------#
# functions (must be contained in working directory)
source("simulation_example_functions.R")

#------------------------------------------------------------------------------#
# setup

# select type of simulation example:
# - "S" for "smooth"
# - "D" for "discontinuous"
# - "P" for "poisson"
# - "NI" for "non-isotonic"
type <- "S"
# type <- Sys.getenv("type")

# sample size for training and validation data (reduce to increase speed)
n <- 2000
# n <- as.integer(Sys.getenv("n"))

nValid <- 5000

#------------------------------------------------------------------------------#
# data generation
nFit <- n
N <- nFit + nValid

id <- 250
# id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) # only on HPC cluster
# id runs from 1 to 500
set.seed(n + id) # for reproducibility of training data

x <- runif(nFit, 0, 10)
y <- switch(type,
  "D" = rgamma(nFit, shape = sqrt(x), scale = pmin(pmax(x, 1), 6))
          + 10 * (x >= 5),
  "P" = rpois(nFit, pmin(pmax(x, 1), 6)),
  "NI" = rgamma(nFit, shape = sqrt(x), scale = pmin(pmax(x, 1), 6))
          - 2 *(x > 7),
  "S" = rgamma(nFit, shape = sqrt(x), scale = pmin(pmax(x, 1), 6)),
)

x <- x + runif(nFit, -1e-6, 1e-6)
y <- y + runif(nFit, -1e-6, 1e-6)

set.seed(id) # for reproducibility of validation data

xValid <- runif(nValid, 0, 10)
yValid <- switch(type,
  "D" = rgamma(
    nValid,
    shape = sqrt(xValid),
    scale = pmin(pmax(xValid, 1), 6)
  ) + 10 * (xValid >= 5),
  "P" = rpois(nValid, pmin(pmax(xValid, 1), 6)),
  "NI" = rgamma(nValid, shape = sqrt(xValid), scale = pmin(pmax(xValid, 1), 6))
          - 2 *(xValid > 7),
  "S" = rgamma(nValid, shape = sqrt(xValid), scale = pmin(pmax(xValid, 1), 6))
)
x <- c(x, xValid)
y <- c(y, yValid)
data <- data.frame(y = y, x = x)

# data frame for distributional regression methods using splines
splineData <- cbind(
  y = y,
  as.data.frame(bs(x, knots = 2 * (1:4), Boundary.knots = c(0, 10)))
)
colnames(splineData) <- c("y", paste0("x", seq_len(ncol(splineData) - 1)))

# training/validation split
newdata <- data[(nFit + 1):N, "x", drop = FALSE]
newSpData <- splineData[(nFit + 1):N, -1]
newY <- y[(nFit + 1):N]
data <- data[1:nFit, ]
splineData <- splineData[1:nFit, ]

#------------------------------------------------------------------------------#
# model fits, prediction, crps
method <- c("idr", "idrbag", "np", "sqr", "tram", "qrf")
crpsRes <- seq_along(method)
names(crpsRes) <- method

dataList <- list(
  # training data for each method
  data,
  data,
  switch(type,"P" = mutate(data, y = factor(y, ordered = TRUE)), data),
  splineData,
  switch(
    type,
    "P" = mutate(splineData, y = factor(y, ordered = TRUE)),
    splineData
  ),
  data
)
MoreArgs <- list(
    # arguments for distributional regression methods
      # idr
      # idrbag
      # np
  bwtype = "adaptive_nn",
  nmulti = 4,
  oykertype = "liracine",
      # tram
      # rq
  tau = seq(0.005, 0.995, 0.001),
      # qrf
  min.node.size = 40,
  sample.fraction = 0.5,
  honesty.fraction = 0.5,
  honesty = TRUE,
  quantiles = seq(0.01, 0.99, 0.01)
)
MoreArgsPred <- list(
  # arguments for predict methods
    digits = 6,                                   # round predictive cdfs
    b = 100,                                      # number of bagging samples
    y = y,                                        # for idr bagging
    K = 5000,                                     # grid for cdfs (tram)
    grd = switch(type,
      "P" = factor(0:max(y), ordered = TRUE),     # grid for np predictions
      seq(min(newY) - 5, max(newY) + 5, 0.1)
    ),
    tau = seq(0.005, 0.995, 0.001),               # for qrf
    range = range(data$x)
)
newdataList <- list(
  # data for out of sample predictions
  newdata,
  newdata,
  newdata,
  newSpData,
  newSpData,
  newdata
)

#------------------------------------------------------------------------------#
# run estimation, prediction, compute crps
for (j in seq_along(method)) {
  fit <- do.call(
    what = probFit,
    args = c(list(method = method[[j]], data = dataList[[j]][1:2000, ]), MoreArgs)
  )
  pred <- do.call(
    what = probPredict,
    c(list(fit = fit, newdata = newdataList[[j]]), MoreArgsPred)
  )
  if (type == "P" & method[j] == "np") {
    pred <- structure(
      map(pred, ~mutate(., points = as.numeric(as.character(points)))),
      class = c("idr", "nppred")
    )
  }
  crpsRes[j] <- mean(crps(pred, newY))
}

# add ideal forecast
crpsRes <- c(
  crpsRes,
  ideal = switch(type,
    "P" = mean(scoringRules::crps_pois(
      y = newY,
      lambda = pmin(pmax(newdata$x, 1), 6)
    )),
    "D" = mean(scoringRules::crps_gamma(
      y = newY - 10 * (newdata$x >= 5),
      shape = sqrt(newdata$x),
      scale = pmin(pmax(newdata$x, 1), 6)
    )),
    "S" = mean(scoringRules::crps_gamma(
      y = newY,
      shape = sqrt(newdata$x),
      scale = pmin(pmax(newdata$x, 1), 6)
    )),
    "NI" = mean(scoringRules::crps_gamma(
      y = newY + 2 *(newdata$x > 7),
      shape = sqrt(newdata$x),
      scale = pmin(pmax(newdata$x, 1), 6)
    ))
  )
)

# print results
crpsRes

# export
save(
  list = "crpsRes",
  file = paste0("iqr_", id, "_", n, "_", type, ".rda")
)
