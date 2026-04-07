.log1mexp <- function(x) {
  stopifnot(all(x <= 0))
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

.log_diff <- function(log_x, log_y) {
  result <- log_x + .log1mexp(log_y - log_x)
  result[is.infinite(log_x) & log_x == -Inf] <- -Inf
  return(result)
}

.log_sum <- function(x, y) {
  ma     <- pmax(x, y)
  result <- ma + log1p(exp(pmin(x, y) - ma))
  result[is.infinite(ma) & ma == -Inf] <- -Inf
  return(result)
}

.check_truncpois_bounds <- function(a, b) {
  if (!is.numeric(a) || any(a < 0L))
    stop("'a' must be non-negative and numeric", call. = FALSE)
  if (!is.numeric(b) || any(b <= 0))
    stop("'b' must be strictly positive and numeric", call. = FALSE)
  if (any(a >= b))
    stop("Lower bound 'a' must be less than upper bound 'b'", call. = FALSE)
  invisible(NULL)
}

.check_lambda <- function(lambda) {
  if (!is.numeric(lambda) || any(lambda <= 0))
    stop("'lambda' must be positive and numeric", call. = FALSE)
  invisible(NULL)
}

.log_denom_truncpois <- function(a, b, lambda) {
  .log_diff(
    stats::ppois(b, lambda, log.p = TRUE, lower.tail = TRUE),
    stats::ppois(a, lambda, log.p = TRUE, lower.tail = TRUE)
  )
}
