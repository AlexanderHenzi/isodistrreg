test_that("idr() single-dimensional not isotone, hazard-rate order", {
  fit <- idr(c(1, 3, 2, 4), data.frame(x = 1:4), stoch = "hazard")
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
