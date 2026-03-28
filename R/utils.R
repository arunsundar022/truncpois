#' Numerically stable log(1 - exp(x))
#'
#' Evaluate \eqn{\log(1 - e^x)} for \eqn{x \le 0} using a two-branch
#' strategy that avoids catastrophic cancellation:
#' \eqn{x > -\log 2} uses \code{log(-expm1(x))}; otherwise
#' \code{log1p(-exp(x))}.
#'
#' @param x Numeric vector with all elements \eqn{\le 0}.
#'
#' @return Numeric vector the same length as \code{x}.
#' @keywords internal
.log1mexp <- function(x) {
  stopifnot(all(x <= 0))
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

#' Numerically stable log-scale subtraction
#'
#' Compute \eqn{\log(x - y)} given \eqn{\log x} and \eqn{\log y} without
#' leaving the log domain. Returns \code{-Inf} wherever \code{log_x} is
#' \code{-Inf}.
#'
#' @param log_x Numeric vector of log-scale values (\eqn{\log x}).
#' @param log_y Numeric vector of log-scale values (\eqn{\log y}), the same
#'   length as \code{log_x}. Must satisfy \code{log_y <= log_x} elementwise.
#'
#' @return Numeric vector the same length as \code{log_x}.
#' @keywords internal
.log_diff <- function(log_x, log_y) {
  result <- log_x + .log1mexp(log_y - log_x)
  result[is.infinite(log_x) & log_x == -Inf] <- -Inf
  return(result)
}

#' Numerically stable log-scale addition
#'
#' Compute \eqn{\log(x + y)} given \eqn{\log x} and \eqn{\log y} via the
#' log-sum-exp trick: \eqn{\max + \log(1 + \exp(\min - \max))}. Returns
#' \code{-Inf} wherever both inputs are \code{-Inf}.
#'
#' @param x Numeric vector of log-scale values.
#' @param y Numeric vector of log-scale values, the same length as \code{x}.
#'
#' @return Numeric vector the same length as \code{x}.
#' @keywords internal
.log_sum <- function(x, y) {
  ma     <- pmax(x, y)
  result <- ma + log1p(exp(pmin(x, y) - ma))
  result[is.infinite(ma) & ma == -Inf] <- -Inf
  return(result)
}

#' Validate truncated Poisson bounds
#'
#' Check that the truncation bounds \code{a} and \code{b} are numeric, satisfy
#' non-negativity / positivity constraints, and that \code{a < b}. Throws an
#' informative error and returns \code{NULL} invisibly on success.
#'
#' @param a Lower truncation bound. Must be a non-negative numeric scalar.
#' @param b Upper truncation bound. Must be a strictly positive numeric scalar
#'   greater than \code{a}.
#'
#' @return \code{NULL} invisibly.
#' @keywords internal
.check_truncpois_bounds <- function(a, b) {
  if (!is.numeric(a) || any(a < 0L))
    stop("'a' must be non-negative and numeric", call. = FALSE)
  if (!is.numeric(b) || any(b <= 0))
    stop("'b' must be strictly positive and numeric", call. = FALSE)
  if (any(a >= b))
    stop("Lower bound 'a' must be less than upper bound 'b'", call. = FALSE)
  invisible(NULL)
}

#' Validate the Poisson rate parameter
#'
#' Check that \code{lambda} is a single, strictly positive numeric scalar.
#' Throws an informative error on failure and returns \code{NULL} invisibly on
#' success.
#'
#' @param lambda Object to validate as the Poisson mean.
#'
#' @return \code{NULL} invisibly.
#' @keywords internal
.check_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0)
    stop("'lambda' must be a single positive numeric", call. = FALSE)
  invisible(NULL)
}

#' Log-scale normalising constant for the truncated Poisson
#'
#' Compute \eqn{\log\!\bigl(P(a < X \le b)\bigr)} for an untruncated
#' \eqn{\text{Poisson}(\lambda)} random variable \eqn{X}, which serves as the
#' log-denominator in all truncated-distribution calculations.
#'
#' @param a Lower truncation bound, already shifted by \code{-1L} relative to
#'   the user-facing \code{a} (i.e. the exclusive lower bound).
#' @param b Upper truncation bound (inclusive), on the original scale.
#' @param lambda Positive numeric scalar, the Poisson mean.
#'
#' @return A single numeric value on the log scale.
#' @keywords internal
.log_denom_truncpois <- function(a, b, lambda) {
  .log_diff(
    stats::ppois(b, lambda, log.p = TRUE, lower.tail = TRUE),
    stats::ppois(a, lambda, log.p = TRUE, lower.tail = TRUE)
  )
}
