test_that("plottruncpois pmf runs without error and returns expected structure", {
  result <- plottruncpois(lambda = 5, a = 2, b = 10, type = "pmf")
  expect_type(result, "list")
  expect_named(result, c("x", "truncated", "untruncated"))
  expect_equal(sum(result$truncated), 1, tolerance = 1e-8)
})

test_that("plottruncpois pmf without comparison omits untruncated values", {
  result <- plottruncpois(lambda = 5, a = 2, b = 10, type = "pmf", compare = FALSE)
  expect_null(result$untruncated)
})

test_that("plottruncpois cdf runs without error and reaches 1", {
  result <- plottruncpois(lambda = 4, a = 1, b = 9, type = "cdf")
  expect_equal(max(result$truncated), 1, tolerance = 1e-8)
})

test_that("plottruncpois quantile type runs without error and is monotonic", {
  result <- plottruncpois(lambda = 5, a = 2, b = 10, type = "quantile")
  expect_true(all(diff(result$truncated) >= 0))
  expect_true(all(result$truncated >= 2 & result$truncated <= 10))
})

test_that("plottruncpois quantile type without comparison omits untruncated values", {
  result <- plottruncpois(lambda = 5, a = 2, b = 10, type = "quantile", compare = FALSE)
  expect_null(result$untruncated)
})
