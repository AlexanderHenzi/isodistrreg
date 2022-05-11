
<!-- README.md is generated from README.Rmd. Please edit that file -->

# isodistrreg: Isotonic Distributional Regression

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/AlexanderHenzi/isodistrreg.svg?token=1kfJSpfJj96s5n1DwsiP&branch=master)](https://travis-ci.com/AlexanderHenzi/isodistrreg)
<!-- badges: end -->

Isotonic distributional regression (IDR) is a powerful nonparametric
technique for the estimation of distributions of a binary or numeric
response variable conditional on numeric or ordinal covariates. IDR
assumes that there is a monotone relationship between the response
variable and the covariates, where the partial order on the covariate
space can be specified by the user, and has no tuning parameters. It can
be used to generate calibrated probabilistic weather forecasts from
ensemble forecasts and observations, and serve as a benchmark in many
other prediction problems.

## Installation

Download the development version from [GitHub](https://github.com/) with

``` r
# install.packages("devtools")
devtools::install_github("AlexanderHenzi/isodistrreg", build_vignettes = TRUE)
```

or install the package from CRAN
(<https://CRAN.R-project.org/package=isodistrreg>):

``` r
install.packages("isodistrreg")
```

If the installation fails on Linux because of the error *Vignette
re-building failed*, try installing with `build_vignettes = FALSE` or
make sure that texinfo is installed (`sudo apt-get install texinfo`).

## Usage Examples

The following basic example illustrates how to use IDR to calibrate
forecasts of accumulated precipitation.

``` r
library(isodistrreg)

# Prepare dataset: Half of the data as training dataset, other half for validation.
# Consult the R documentation (?rain) for details about the dataset.
data(rain)
trainingData <- subset(rain, date <= "2012-01-09")
validationData <- subset(rain, date > "2012-01-09")

# Variable selection: use HRES and the perturbed forecasts P1, ..., P50
varNames <- c("HRES", paste0("P", 1:50))

# Partial orders on variable groups: Usual order of numbers on HRES (group '1') and
# increasing convex order on the remaining variables (group '2').
groups <- setNames(c(1, rep(2, 50)), varNames)
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

Henzi, A., Ziegel, J.F. and Gneiting, T. (2021), Isotonic distributional
regression. J R Stat Soc Series B, 83: 963-993.
<https://doi.org/10.1111/rssb.12450>

Henzi, A., Moesching, A. & Duembgen, L. Accelerating the
Pool-Adjacent-Violators Algorithm for Isotonic Distributional
Regression. Methodol Comput Appl Probab (2022).
<https://doi.org/10.1007/s11009-022-09937-2>

The dataframe `precipData_caseStudy.rda` contains all data for the case
study in the paper.
