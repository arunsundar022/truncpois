test_that("functions behave correctly for extreme lambda", {
  lambda <- 1e5
  a <- lambda - 500
  b <- lambda + 500

  d <- dtruncpois(a:b, lambda, a = a, b = b, log = TRUE)
  expect_false(any(is.nan(d)))
  expect_true(all(is.finite(d)))

  p <- ptruncpois(seq(a, b, by = 100), lambda, a = a, b = b)
  expect_false(any(is.nan(p)))
  expect_true(all(diff(p) >= 0))

  probs <- c(0.1, 0.5, 0.9)
  q <- qtruncpois(probs, lambda, a = a, b = b)
  cdf_at_q <- ptruncpois(q, lambda, a = a, b = b)
  expect_true(all(cdf_at_q >= probs - 1e-8))

  set.seed(11)
  x <- rtruncpois(500, lambda, a = a, b = b, method = "bounded")
  expect_true(all(x >= a & x <= b))
})

test_that("functions behave correctly for the narrowest possible truncation window", {
  lambda <- 5
  a <- 3
  b <- a + 1L 

  probs <- dtruncpois(c(a, b), lambda, a = a, b = b)
  expect_equal(sum(probs), 1, tolerance = 1e-10)
  expect_equal(dtruncpois(a - 1L, lambda, a = a, b = b), 0)
  expect_equal(dtruncpois(b + 1L, lambda, a = a, b = b), 0)

  for (method in c("direct", "inversion", "bounded")) {
    set.seed(21)
    x <- rtruncpois(200, lambda, a = a, b = b, method = method)
    expect_true(all(x %in% c(a, b)), info = paste("method:", method))
  }
})

test_that("rtruncpois large-sample generation converges across all methods", {
  lambda <- 4
  a <- 1
  b <- 10
  n <- 1e5

  for (method in c("direct", "inversion", "bounded")) {
    set.seed(31)
    x <- rtruncpois(n, lambda, a = a, b = b, method = method)
    expect_true(all(x >= a & x <= b), info = paste("method:", method))
    expect_equal(mean(x), extruncpois(lambda, a = a, b = b), tolerance = 0.02,
                 label = paste("method:", method))
    expect_equal(var(x), vartruncpois(lambda, a = a, b = b), tolerance = 0.05,
                 label = paste("method:", method))
  }
})

test_that("rtruncpois respects bounds across many random seeds", {
  lambda <- 3
  a <- 0
  b <- 8

  for (seed in 1:20) {
    set.seed(seed)
    for (method in c("direct", "inversion", "bounded")) {
      x <- rtruncpois(100, lambda, a = a, b = b, method = method)
      expect_true(all(x >= a & x <= b),
                  info = paste("seed:", seed, "method:", method))
    }
  }
})
