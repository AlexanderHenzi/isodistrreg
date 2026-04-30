test_that("idr() single-dimensional not isotone", {
  fit <- idr(c(1, 3, 2, 4), data.frame(x = 1:4))
  preds <- predict(fit, data.frame(x = 1:4))
  expect_equal(
    preds$cdf,
    t(
      matrix(
        c(
          1, 1, 1, 1,
          0, 0.5, 1, 1,
          0, 0.5, 1, 1,
          0, 0, 0, 1
        ),
        nrow = 4,
        ncol = 4
      )
    )
  )
})

test_that("idr() single-dimensional not isotone 2", {
  fit <- idr(c(3, 5, 3), data.frame(x = c(7, 11, 13)))
  preds <- predict(fit, data.frame(x = c(7, 11, 13)))
  expect_equal(
    preds$cdf,
    matrix(
      c(
        1, 0.5, 0.5,
        1, 1, 1
      ),
      nrow = 3,
      ncol = 2
    )
  )
})

test_that("idr() single-dimensional single covariate", {
  fit <- idr(c(3, 5), data.frame(x = rep(7, 2)), weights = c(1, 2))
  preds <- predict(fit, data.frame(x = 7))
  expect_equal(
    preds$cdf,
    matrix(
      c(1 / 3, 1),
      nrow = 1,
      ncol = 2
    )
  )
})

test_that("idr() single-dimensional single response", {
  fit <- idr(c(3, 3), data.frame(x = c(5, 7)), weights = c(1, 2))
  preds <- predict(fit, data.frame(x = c(5, 7)))
  expect_equal(
    preds$cdf,
    matrix(
      c(1, 1),
      nrow = 2,
      ncol = 1
    )
  )
})
