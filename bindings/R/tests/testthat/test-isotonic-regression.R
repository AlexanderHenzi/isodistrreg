test_that("isotonic_regression with covariates returns isotone means", {
  result <- isotonic_regression(c(2, 3, 1, 4, 5), X = as.double(1:5))
  # Result should be non-decreasing
  expect_true(all(diff(result) >= 0))
  # Already isotone except position 3 (value 1), which gets pooled with
  # neighbors
  expect_equal(result, c(2, 2, 2, 4, 5))
})

test_that("isotonic_regression already isotone input is unchanged", {
  result <- isotonic_regression(c(1, 2, 3, 4), X = as.double(1:4))
  expect_equal(result, c(1, 2, 3, 4))
})

test_that("isotonic_regression with weights", {
  result <- isotonic_regression(
    c(3, 2, 4, 1),
    X = as.double(1:4),
    weights = c(1, 2, 1, 1)
  )
  expect_true(all(diff(result) >= 0))
})

test_that("isotonic_regression pre-sorted (no covariates)", {
  result <- isotonic_regression(sort(c(3, 1, 2, 5)))
  expect_equal(result, c(1, 2, 3, 5))
})

test_that("isotonic_regression pre-sorted with weights", {
  result <- isotonic_regression(sort(c(2, 1, 3)), weights = c(1, 2, 1))
  expect_true(all(diff(result) >= 0))
})

test_that("isotonic_regression constant input stays constant", {
  result <- isotonic_regression(rep(5, 4), X = as.double(1:4))
  expect_equal(result, rep(5, 4))
})

test_that("isotonic_regression decreasing direction", {
  result <- isotonic_regression(
    c(5, 3, 4, 1),
    X = as.double(1:4),
    decreasing = TRUE
  )
  # Result should be non-increasing
  expect_true(all(diff(result) <= 0))
})

test_that("isotonic_regression single element", {
  result <- isotonic_regression(c(42), X = c(1))
  expect_equal(result, 42)
})
