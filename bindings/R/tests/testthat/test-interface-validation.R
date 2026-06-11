test_that("crps works for fits with few jump points and recycles scalar y", {
  fit <- idr(c(0, 0, 1, 1), data.frame(x = rep(1, 4)), progress = FALSE)
  pred <- predict(fit, data.frame(x = 1:2))
  # Two jump points, two predictions: used to fail with "subscript out of bounds".
  expect_equal(crps(pred, c(0, 1)), c(0.25, 0.25))

  fit5 <- idr(as.double(1:5), data.frame(x = as.double(1:5)), progress = FALSE)
  pred5 <- predict(fit5, data.frame(x = 1:5))
  # Scalar y must be recycled to all predictions, not score only the first.
  expect_equal(crps(pred5, 3), c(2, 1, 0, 1, 2))

  fit1 <- idr(rep(1, 3), data.frame(x = as.double(1:3)), progress = FALSE)
  pred1 <- predict(fit1, data.frame(x = 1:3))
  expect_length(crps(pred1, c(1, 1, 1)), 3)
})

test_that("ordered factor covariates are converted by level, not by label", {
  X <- data.frame(f = factor(c("lo", "mid", "hi", "hi"),
    levels = c("lo", "mid", "hi"), ordered = TRUE
  ))
  fit <- idr(c(1, 2, 3, 3), X, progress = FALSE)
  expect_s3_class(fit, "idrfit")
  pred <- predict(fit, X)
  expect_equal(unname(pred$cdf), unname(fit$cdf))
})

test_that("prediction inputs are validated", {
  fit <- idr(c(0, 0, 1, 1), data.frame(x = rep(1, 4)), progress = FALSE)
  pred <- predict(fit, data.frame(x = 1))
  expect_error(cdf(pred, NA_real_), "must not contain NAs")
  expect_error(qpred(pred, NA_real_), "no NAs")
  # Integer vectors are documented as numeric and must be coerced.
  expect_silent(cdf(pred, 0L))
  expect_silent(qpred(pred, c(0L, 1L)))
})

test_that("invalid fit inputs give targeted errors instead of Rust panics", {
  X <- data.frame(x = as.double(1:4))
  expect_error(idr(c(1, 2, Inf, 4), X, progress = FALSE), "finite")
  expect_error(
    idr(c(1, 2, 3, 4), X, weights = c(1, NA, 1, 1), progress = FALSE),
    "finite non-negative"
  )
  expect_error(
    idr(c(1, 2, 3, 4), X, weights = rep(0, 4), progress = FALSE),
    "at least one weight"
  )
  expect_error(
    idr(c(1, 2, 3, 4), X, y_observed = rep(0, 4), progress = FALSE),
    "uncensored"
  )
  expect_error(
    idr(c(1, 2, 3, 4), X, y_observed = c(2, 1, 1, 1), progress = FALSE),
    "TRUE/1 or FALSE/0"
  )
  expect_error(
    idr(c(1, 2, 3, 4), X, pars = list(bogus = 1), progress = FALSE),
    "bogus"
  )
  X2 <- data.frame(a = as.double(1:4), b = as.double(4:1))
  expect_error(
    idr(c(1, 2, 3, 4), X2,
      groups = setNames(c(1, 1, 1), c("a", "b", "zzz")),
      orders = c(comp = 1), progress = FALSE
    ),
    "not in 'X'"
  )
  expect_error(
    idr(c(1, 2, 3, 4), X, orders = c(1), progress = FALSE),
    "named"
  )
})

test_that("partial pars lists are merged with the documented defaults", {
  X2 <- data.frame(a = as.double(1:4), b = as.double(4:1))
  fit <- idr(c(2, 1, 3, 5), X2, pars = list(verbose = FALSE), progress = FALSE)
  expect_s3_class(fit, "idrfit")
  fit_null <- idr(c(2, 1, 3, 5), X2, pars = NULL, progress = FALSE)
  expect_s3_class(fit_null, "idrfit")
})

test_that("predict matches columns by name and supports in-sample prediction", {
  X2 <- data.frame(a = as.double(1:5), b = c(0, 0, 10, 10, 10))
  fit <- idr(as.double(1:5), X2, progress = FALSE)
  p1 <- predict(fit, data.frame(a = 3, b = 10))
  p2 <- predict(fit, data.frame(b = 10, a = 3))
  expect_identical(p1$cdf, p2$cdf)
  expect_error(predict(fit, data.frame(foo = 3, bar = 10)), "same variables")
  # In-sample predictions when data is omitted, as documented.
  p_in <- predict(fit)
  expect_equal(unname(p_in$cdf), unname(fit$cdf))
})

test_that("idrbag returns predictions of class idr, as documented", {
  set.seed(7)
  y <- rnorm(20)
  X <- data.frame(x = rnorm(20))
  nd <- data.frame(x = rnorm(5))
  bag <- idrbag(y, X, newdata = nd, b = 2, p = 0.5, progress = FALSE)
  expect_s3_class(bag, "idr")
  expect_equal(nrow(cdf(bag, 0)), 5)
  expect_length(crps(bag, rnorm(5)), 5)

  grid <- c(-1, 0, 1)
  bag_grid <- idrbag(y, X,
    newdata = nd, b = 2, p = 0.5, grid = grid,
    progress = FALSE
  )
  expect_identical(bag_grid$points, grid)
  expect_identical(dim(bag_grid$cdf), c(5L, 3L))

  expect_error(
    idrbag(y, X, newdata = nd, b = 2, p = 0.5, grid = c(1, NA), progress = FALSE),
    "without NAs"
  )
  expect_error(
    idrbag(y, X, newdata = nd, b = NA_real_, p = 0.5, progress = FALSE),
    "positive integer"
  )
})

test_that("pit is randomized only at jump points of the predictive CDF", {
  fit <- idr(c(0, 0, 1, 1), data.frame(x = rep(1, 4)), progress = FALSE)
  pred <- predict(fit, data.frame(x = c(1, 1)))
  # y strictly between jumps: the CDF is continuous there, PIT is deterministic.
  set.seed(99)
  expect_equal(pit(pred, c(1.3, 0.3)), c(1.0, 0.5))
  # y at a jump: the PIT must be randomized within (F(y-), F(y)].
  set.seed(1)
  vals <- pit(pred, c(1, 1))
  expect_true(all(vals >= 0.5 & vals <= 1))
})

test_that("isotonic_regression validates input and accepts integer vectors", {
  expect_equal(isotonic_regression(c(3L, 1L, 2L)), c(2, 2, 2))
  expect_error(isotonic_regression(c(1, NA, 3)), "finite")
  expect_error(isotonic_regression(c(3, 2, 1), weights = c(1, -1, 1)), "non-negative")
  expect_error(isotonic_regression(c(1, 2, 3), X = c(1, 2)), "equal length")
  expect_error(isotonic_regression(c(3, 2, 1), weights = c(0, 0, 0)), "positive")
})
