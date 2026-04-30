test_that("cdf.data.frame computes empirical CDF at thresholds", {
  predictions <- data.frame(a = c(1, 3), b = c(2, 4), c = c(3, 5))
  result <- cdf(predictions, thresholds = c(0, 2, 3, 6))

  # Row 1: values {1, 2, 3}. ecdf: P(X <= 0) = 0, P(X <= 2) = 2/3,
  #   P(X <= 3) = 1, P(X <= 6) = 1
  # Row 2: values {3, 4, 5}. ecdf: P(X <= 0) = 0, P(X <= 2) = 0,
  #   P(X <= 3) = 1/3, P(X <= 6) = 1
  expected <- matrix(
    c(
      0, 2 / 3, 1, 1,
      0, 0, 1 / 3, 1
    ),
    nrow = 2,
    ncol = 4,
    byrow = TRUE
  )
  expect_equal(result, expected)
})

test_that("qpred.data.frame computes empirical quantiles", {
  predictions <- data.frame(a = c(1, 10), b = c(2, 20), c = c(3, 30))
  result <- qpred(predictions, quantiles = c(0, 0.5, 1))

  # stats::quantile type 1 (inverse of ecdf):
  # Row 1: values {1, 2, 3}. Q(0) = 1, Q(0.5) = 2, Q(1) = 3
  # Row 2: values {10, 20, 30}. Q(0) = 10, Q(0.5) = 20, Q(1) = 30
  expected <- matrix(
    c(
      1, 2, 3,
      10, 20, 30
    ),
    nrow = 2,
    ncol = 3,
    byrow = TRUE
  )
  expect_equal(result, expected)
})
