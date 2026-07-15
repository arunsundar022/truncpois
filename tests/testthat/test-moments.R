test_that("extruncpois zero-truncated mean is greater than lambda", {
  lambda <- 2
  expect_gt(extruncpois(lambda, a = 1), lambda)
})

test_that("extruncpois with default bounds equals lambda", {
  expect_equal(extruncpois(5), 5, tolerance = 1e-10)
  expect_equal(extruncpois(1), 1, tolerance = 1e-10)
})

test_that("extruncpois converges to lambda as b increases", {
  means <- sapply(5:30, function(b) extruncpois(lambda = 3, a = 0, b = b))
  expect_equal(means[length(means)], 3, tolerance = 1e-4)
})

test_that("extruncpois stops for invalid lambda", {
  expect_error(extruncpois(-1))
  expect_error(extruncpois(0))
})

test_that("vartruncpois is non-negative", {
  expect_gte(vartruncpois(4, a = 0, b = 8), 0)
  expect_gte(vartruncpois(1, a = 1, b = 5), 0)
})

test_that("vartruncpois with default bounds equals lambda (Poisson variance = mean)", {
  expect_equal(vartruncpois(5), 5, tolerance = 1e-8)
})

test_that("truncation reduces variance", {
  lambda <- 5
  untrunc_var <- lambda
  trunc_var <- vartruncpois(lambda, a = 2, b = 8)
  expect_lt(trunc_var, untrunc_var)
})

test_that("vartruncpois large-sample convergence", {
  set.seed(7)
  s <- rtruncpois(10000, lambda = 4, a = 1, b = 10)
  expect_equal(var(s), vartruncpois(4, a = 1, b = 10), tolerance = 0.1)
})

test_that("medtruncpois result is within [a, b]", {
  m <- medtruncpois(lambda = 5, a = 2, b = 9)
  expect_gte(m, 2)
  expect_lte(m, 9)
})

test_that("medtruncpois equals qtruncpois at p = 0.5", {
  lambda <- 4
  a <- 1
  b <- 12
  expect_equal(
    medtruncpois(lambda, a, b),
    qtruncpois(0.5, lambda, a, b)
  )
})

test_that("medtruncpois with default bounds is consistent with qpois", {
  expect_equal(
    medtruncpois(lambda = 3),
    qtruncpois(0.5, lambda = 3)
  )
})

test_that("modtruncpois issues warning and returns two values for integer lambda", {
  expect_warning(modtruncpois(lambda = 3, a = 0, b = 10), "mode is not unique")
  expect_warning(modtruncpois(lambda = 5, a = 0, b = 10), "mode is not unique")
  suppressWarnings({
    expect_equal(modtruncpois(lambda = 3, a = 0, b = 10), c(2L, 3L))
    expect_equal(modtruncpois(lambda = 5, a = 0, b = 10), c(4L, 5L))
  })
})

test_that("modtruncpois returns unique mode for non-integer lambda", {
  m <- modtruncpois(lambda = 2.5, a = 0, b = 10)
  expect_length(m, 1)
  expect_equal(m, 2L)
})

test_that("modtruncpois clamps result to [a, b]", {
  m <- modtruncpois(lambda = 1, a = 3, b = 10)
  expect_equal(m, 3L)
})

test_that("modtruncpois stops for invalid inputs", {
  expect_error(modtruncpois(lambda = -1))
  expect_error(modtruncpois(lambda = 3, a = 5, b = 3))
})
