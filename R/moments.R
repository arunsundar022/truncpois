#' Truncated Poisson mean
#'
#' Compute the expected value of a truncated Poisson distribution,
#' \eqn{E[X \mid a < X \le b]}, evaluated in log-space for numerical
#' stability.
#'
#' @param lambda Positive numeric scalar, the Poisson mean.
#' @param a Lower truncation bound (inclusive). Default is \code{0}.
#' @param b Upper truncation bound (inclusive). Default is \code{Inf}.
#'
#' @return A single positive numeric value.
#' @export
extruncpois <- function(lambda, a = 0L, b = Inf) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)

  a         <- a - 1L
  log_denom <- .log_denom_truncpois(a, b, lambda)
  log_num   <- .log_diff(
    stats::ppois(b - 1L, lambda, log.p = TRUE, lower.tail = TRUE),
    stats::ppois(a - 1L, lambda, log.p = TRUE, lower.tail = TRUE)
  )
  return(exp(log(lambda) + log_num - log_denom))
}

#' Truncated Poisson variance
#'
#' Compute the variance of a truncated Poisson distribution,
#' \eqn{\text{Var}(X \mid a < X \le b)}, via the identity
#' \eqn{E[X^2] - E[X]^2}. The second moment is obtained by expressing
#' \eqn{E[X^2] = \lambda^2 \cdot P(a{+}2 < X \le b{+}2) / P(a < X \le b) + E[X]},
#' keeping all CDF evaluations on the log scale.
#'
#' @param lambda Positive numeric scalar, the Poisson mean.
#' @param a Lower truncation bound (inclusive). Default is \code{0}.
#' @param b Upper truncation bound (inclusive). Default is \code{Inf}.
#'
#' @return A single non-negative numeric value.
#'
#' @seealso \code{\link{extruncpois}} for the corresponding mean.
#' @export
vartruncpois <- function(lambda, a = 0L, b = Inf) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)

  a           <- a - 1L
  log_denom   <- .log_denom_truncpois(a, b, lambda)
  mu          <- extruncpois(lambda, a + 1L, b)
  log_fac_num <- .log_diff(
    stats::ppois(b - 2L, lambda, log.p = TRUE, lower.tail = TRUE),
    stats::ppois(a - 2L, lambda, log.p = TRUE, lower.tail = TRUE)
  )

  ex2 <- exp(2 * log(lambda) + log_fac_num - log_denom) + mu
  return(ex2 - mu^2)
}
