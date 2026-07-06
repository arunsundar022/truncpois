test_that("mletruncpois recovers a known lambda from a large sample", {
  set.seed(1)
  x <- rtruncpois(5000, lambda = 6, a = 1)
  fit <- mletruncpois(x, a = 1)
  expect_equal(fit$lambda, 6, tolerance = 0.1)
  expect_equal(fit$n, 5000)
})

test_that("mletruncpois recovers lambda for a doubly-truncated sample", {
  set.seed(2)
  x <- rtruncpois(5000, lambda = 8, a = 2, b = 15)
  fit <- mletruncpois(x, a = 2, b = 15)
  expect_equal(fit$lambda, 8, tolerance = 0.15)
})

test_that("mletruncpois stops for non-integer x", {
  expect_error(mletruncpois(c(1.5, 2, 3), a = 1), "must contain only integers")
})

test_that("mletruncpois stops for x outside [a, b]", {
  expect_error(mletruncpois(c(0, 1, 2), a = 1), "must lie within")
})

test_that("mletruncpois stops for empty x", {
  expect_error(mletruncpois(numeric(0)), "non-empty")
})
