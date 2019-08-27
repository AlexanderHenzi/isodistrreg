
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isodistrreg: Isotonic Distributional Regression

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/AlexanderHenzi/isodistrreg.svg?branch=master)](https://travis-ci.org/AlexanderHenzi/isodistrreg)
<!-- badges: end -->

## Authors

The package is based on the paper XXX and written and maintained by
Alexander Henzi.

## Installation

Download the development version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("AlexanderHenzi/isodistrreg")
```

## Usage Examples

The following basic example shows how to use IDR for precipitation
forecast calibration.

``` r
library(isodistrreg)

# Prepare dataset: Half of the data as training dataset, other half for validation.
# Consult the R documentation (?precip) for details about the dataset.
data(precip)
trainingData <- subset(precip, dates <= "2012-01-09")
validationData <- subset(precip, dates > "2012-01-09")

# Variable selection: Use all forecasts (high resolution HRES, perturbed
# forecasts P1, ..., P50 and control run CNT for the perturbed forecasts)
varNames <- c("HRES", "CNT", paste0("P", 1:50))

# Partial orders on variable groups: Usual order of numbers on HRES (group '1') and
# increasing convex order on the remaining variables (group '2').
# See the reference for more details on how to select an appropriate partial order.
groups <- setNames(c(1, rep(2, 51)), varNames)
orders <- c("comp" = 1, "icx" = 2)

# Fit IDR to training dataset.
fit <- idr(
  y = trainingData[["obs"]],
  X = trainingData[, varNames],
  groups = groups,
  orders = orders
)

# Make prediction for the first day in the validation data:
firstPrediction <- predict(fit, data = validationData[1, varNames])
plot(firstPrediction)

# Use cdf() and qpred() to make probability and quantile forecasts:

## What is the probability of precipitation?
1 - cdf(firstPrediction, thresholds = 0)

## What are the predicted 10%, 50% and 90% quantiles for precipitation?
qpred(firstPrediction, quantiles = c(0.1, 0.5, 0.9))

# Make predictions for the complete verification dataset and compare IDR calibrated
# forecasts to the raw ensemble (ENS):
predictions <- predict(fit, data = validationData[, varNames])
y <- validationData[["obs"]]

## Continuous ranked probability score (CRPS):
CRPS <- cbind(
  "ens" = crps(validationData[, varNames], y),
  "IDR" = crps(predictions, y)
)
apply(CRPS, 2, mean)

## Brier score for probability of precipitation:
BS <- cbind(
  "ens" = bscore(validationData[, varNames], thresholds = 0, y),
  "IDR" = bscore(predictions, thresholds = 0, y)
)
apply(BS, 2, mean)

## Quantile score of forecast for 90% quantile:
QS90 <- cbind(
  "ens" = qscore(validationData[, varNames], quantiles = 0.9, y),
  "IDR" = qscore(predictions, quantiles = 0.9, y)
)
apply(QS90, 2, mean)

## Check calibration using (randomized) PIT histograms:
pitEns <- pit(validationData[, varNames], y)
pitIdr <- pit(predictions, y)

par(mfrow = c(1, 2))
hist(pitEns, main = "PIT of raw ensemble forecasts", freq = FALSE)
hist(pitIdr, main = "PIT of IDR calibrated forecasts", freq = FALSE)
```

## References
