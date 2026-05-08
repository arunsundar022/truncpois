test_that("dtruncpois sums to 1 over support", {
  lambda <- 3; a <- 1; b <- 10
  x <- a:b
  expect_equal(sum(dtruncpois(x, lambda, a, b)), 1)
})

test_that("dtruncpois returns 0 outside [a, b]", {
  expect_equal(dtruncpois(0,  lambda = 3, a = 1, b = 10), 0)
  expect_equal(dtruncpois(11, lambda = 3, a = 1, b = 10), 0)
})

test_that("dtruncpois log=TRUE is consistent with log=FALSE", {
  x <- 1:8
  d  <- dtruncpois(x, lambda = 4, a = 1, b = 8, log = FALSE)
  dl <- dtruncpois(x, lambda = 4, a = 1, b = 8, log = TRUE)
  expect_equal(log(d), dl)
})

test_that("dtruncpois with default bounds matches dpois", {
  x <- 0:15
  expect_equal(dtruncpois(x, lambda = 5), stats::dpois(x, lambda = 5),
               tolerance = 1e-12)
})

test_that("dtruncpois issues warning when lambda length mismatches x", {
  expect_warning(dtruncpois(1:4, lambda = c(1, 2, 3), a = 0, b = Inf),
                 "'lambda' length")
})

test_that("dtruncpois stops for non-numeric x", {
  expect_error(dtruncpois("a", lambda = 2), "'x' must be numeric")
})

test_that("ptruncpois reaches exactly 1 at upper bound", {
  expect_equal(ptruncpois(8, lambda = 3, a = 1, b = 8), 1)
})

test_that("ptruncpois is 0 strictly below lower bound", {
  expect_equal(ptruncpois(0, lambda = 3, a = 1, b = 8), 0)
})

test_that("ptruncpois is non-decreasing", {
  q <- 1:8
  cdf <- ptruncpois(q, lambda = 3, a = 1, b = 8)
  expect_true(all(diff(cdf) >= 0))
})

test_that("ptruncpois with default bounds matches ppois", {
  q <- 0:12
  expect_equal(ptruncpois(q, lambda = 5), stats::ppois(q, lambda = 5),
               tolerance = 1e-12)
})

test_that("ptruncpois lower.tail=FALSE is complement of lower.tail=TRUE", {
  q <- 3
  p_lo <- ptruncpois(q, lambda = 4, a = 1, b = 9, lower.tail = TRUE)
  p_hi <- ptruncpois(q, lambda = 4, a = 1, b = 9, lower.tail = FALSE)
  expect_equal(p_lo + p_hi, 1)
})

test_that("ptruncpois log.p=TRUE is consistent with log.p=FALSE", {
  q <- 2:7
  p  <- ptruncpois(q, lambda = 4, a = 1, b = 8, log.p = FALSE)
  pl <- ptruncpois(q, lambda = 4, a = 1, b = 8, log.p = TRUE)
  expect_equal(log(p), pl, tolerance = 1e-12)
})

test_that("qtruncpois is right inverse of ptruncpois", {
  p <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q <- qtruncpois(p, lambda = 5, a = 1, b = 12)
  cdf_at_q <- ptruncpois(q, lambda = 5, a = 1, b = 12)
  expect_true(all(cdf_at_q >= p - 1e-10))
})

test_that("qtruncpois with log.p=TRUE gives same result as log.p=FALSE", {
  p <- c(0.1, 0.5, 0.9)
  q1 <- qtruncpois(p,        lambda = 4, a = 1, b = 10, log.p = FALSE)
  q2 <- qtruncpois(log(p),   lambda = 4, a = 1, b = 10, log.p = TRUE)
  expect_equal(q1, q2)
})

test_that("qtruncpois with default bounds matches qpois", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qtruncpois(p, lambda = 5), stats::qpois(p, lambda = 5))
})

test_that("qtruncpois stops for out-of-range p", {
  expect_error(qtruncpois(-0.1, lambda = 3, a = 0, b = Inf))
  expect_error(qtruncpois(1.1,  lambda = 3, a = 0, b = Inf))
})

test_that("rtruncpois returns values within [a, b]", {
  set.seed(1)
  for (method in c("direct", "inversion", "bounded")) {
    x <- rtruncpois(200, lambda = 4, a = 2, b = 8, method = method)
    expect_true(all(x >= 2 & x <= 8), info = paste("method:", method))
  }
})

test_that("rtruncpois large-sample mean converges to extruncpois", {
  set.seed(42)
  x <- rtruncpois(5000, lambda = 3, a = 1, b = Inf)
  expect_equal(mean(x), extruncpois(3, a = 1), tolerance = 0.05)
})

test_that("rtruncpois returns correct length", {
  expect_length(rtruncpois(10, lambda = 2), 10)
  expect_length(rtruncpois(0,  lambda = 2), 0)
})

test_that("rtruncpois stops for non-scalar n", {
  expect_error(rtruncpois(c(5, 5), lambda = 2))
})
