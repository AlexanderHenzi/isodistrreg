test_that("idrbag() single-dimensional", {
  x <- idrbag(
    c(1, 3, 2, 5), data.frame(x = 1:4),
    b = 100,
    p = 0.5,
    newdata = data.frame(x = c(1, 2.5, 4)),
    progress = FALSE
  )
  expect_equal(length(x$points), 4)
})

test_that("idrbag() single-dimensional censored", {
  x <- idrbag(
    c(1, 3, 2, 4), data.frame(x = 1:4), c(0, 0, 1, 0),
    b = 100,
    p = 0.5,
    newdata = data.frame(x = c(1, 2.5, 4)),
    grid = c(2, 3),
    progress = FALSE
  )
  expect_equal(length(x$points), 2)
})

test_that("idrbag() seed produces identical results across calls", {
  args <- list(
    y        = c(1, 3, 2, 5, 4, 2, 3),
    X        = data.frame(x = c(1, 4, 2, 6, 3, 2, 5)),
    b        = 20,
    p        = 0.6,
    newdata  = data.frame(x = c(1, 3, 5)),
    progress = FALSE,
    seed     = 42
  )
  fit1 <- do.call(idrbag, args)
  fit2 <- do.call(idrbag, args)
  expect_identical(fit1$cdf, fit2$cdf)
})
