#' Truncated Poisson Mean
#'
#' Compute the expected value of a truncated poisson distribution
#'
#' Compute the expected value of a truncated Poisson distribution,
#' \eqn{E[X \mid a \le X \le b]}, evaluated in log-space for numerical
#' stability. With \code{a = 0} and \code{b = Inf} as defaults, the function returns
#' the mean of a non-truncated Poisson distribution.
#'
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Upper truncation bound (inclusive). Default is \code{b = Inf}.
#'
#' @return A single positive numeric value representing the expected value.
#'
#' @note
#' This implementation adapts the continuous truncated distribution framework of
#' Nadarajah and Kotz (2006) to the discrete Poisson setting.
#' Key differences include:
#' \itemize{
#' \item All internal calculations are performed on the log-scale for numerical stability
#' \item The normalizing constant is computed as log-differnce of the cumulative probabilities
#'    rather than densities
#' }
#'
#' @references
#' Nadarajah, S. and Kotz, S. (2006).
#' \dQuote{R Programs for Computing Truncated Distributions}
#' \emph{Journal of Statistical Software}
#' \url{https://www.jstatsoft.org/v16/c02/}
#'
#' @author Arun Sundar Paulraj and Keefe Murphy
#'
#' @export
#'
#' @examples
#' # The mean of a zero-truncated Poisson distribution is strictly greater than lambda
#' lambda_val <- 2
#' extruncpois(lambda_val, a = 1) # > 2
#'
#' # Examine how the mean changes as the upper bound increases
#' bounds <- 1:10
#' means <- sapply(bounds, function(b) extruncpois(lambda = 5, a = 0, b = b))
#'
#' # Plot demonstrating convergence to the untruncated mean
#' plot(bounds, means, type = "b", pch = 19, col = "blue",
#'      ylab = "Expected Value", xlab = "Upper Bound (b)",
#'      main = "Expected Value of Right-Truncated Poisson(lambda=5)")
#' abline(h = 5, lty = 2, col = "red") # Untruncated mean
extruncpois <- function(lambda, a = 0L, b = Inf) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)

  # if (a < 1L) {
  #   qlo <- stats::qpois(.Machine$double.eps, lambda, lower.tail = TRUE,  log.p = FALSE)
  #   qhi <- stats::qpois(.Machine$double.eps, lambda, lower.tail = FALSE, log.p = FALSE)
  #   support <- seq(max(floor(a), qlo), min(floor(b), qhi))
  #   probs <- dtruncpois(support, lambda, a, b, log = TRUE)
  #   w <- exp(probs - max(probs))
  #   return(sum(support * w) / sum(w))
  # } else {
    a_adj <- a - 1L
    log_denom <- .log_denom_truncpois(a_adj, b, lambda)
    log_num <- .log_diff(
      stats::ppois(b - 1L, lambda, log.p = TRUE, lower.tail = TRUE),
      stats::ppois(a_adj - 1L, lambda, log.p = TRUE, lower.tail = TRUE)
    )
    return(exp(log(lambda) + log_num - log_denom))
  #}
}

#' Truncated Poisson Variance
#'
#' Compute the variance of a truncated poisson distribution
#'
#' Compute the variance of a truncated Poisson distribution,
#' \eqn{\text{Var}(X \mid a \le X \le b)}, via the identity
#' \eqn{E[X^2] - E[X]^2}. The second moment is obtained by expressing
#' \eqn{E[X^2] = \lambda^2 \cdot P(a{-}2 \le X \le b{-}2) / P(a \le X \le b) + E[X]},
#' keeping all CDF evaluations on the log scale. With \code{a = 0} and
#' \code{b = Inf} as defaults, the function returns the variance of a
#' non-truncated Poisson distribution.
#'
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Upper truncation bound (inclusive). Default is \code{b = Inf}.
#'
#' @return A single non-negative numeric value representing the variance.
#'
#' @note
#' This implementation adapts the continuous truncated distribution framework of
#' Nadarajah and Kotz (2006) to the discrete Poisson setting.
#' Key differences include:
#' \itemize{
#' \item All internal calculations are performed on the log-scale for numerical stability
#' \item The normalizing constant is computed as log-differnce of the cumulative probabilities
#'    rather than densities
#' }
#'
#' @seealso \code{\link{extruncpois}} for the corresponding mean.
#'
#' @references
#' Nadarajah, S. and Kotz, S. (2006).
#' \dQuote{R Programs for Computing Truncated Distributions}
#' \emph{Journal of Statistical Software}
#' \url{https://www.jstatsoft.org/v16/c02/}
#'
#' @author Arun Sundar Paulraj and Keefe Murphy
#'
#' @export
#'
#' @examples
#' # Poisson distribution has variance == mean.
#' # Truncation reduces the variance (underdispersion).
#' lambda_val <- 5
#' untrunc_var <- lambda_val
#' trunc_var <- vartruncpois(lambda_val, b = 6)
#'
#' cat("Untruncated Variance:", untrunc_var, "\n")
#' cat("Truncated Variance  :", trunc_var, "\n")
#'
#' # Verify computationally with a large sample simulation
#' set.seed(42)
#' samples <- rtruncpois(10000, lambda = 5, b = 6)
#' var(samples) # Should be close to trunc_var
vartruncpois <- function(lambda, a = 0L, b = Inf) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)

  # if (a < 2L) {
  #   qlo <- stats::qpois(.Machine$double.eps, lambda, lower.tail = TRUE,  log.p = FALSE)
  #   qhi <- stats::qpois(.Machine$double.eps, lambda, lower.tail = FALSE, log.p = FALSE)
  #   support <- seq(max(floor(a), qlo), min(floor(b), qhi))
  #   probs <- dtruncpois(support, lambda, a, b, log = TRUE)
  #   w <- exp(probs - max(probs))
  #   mu <- sum(support * w) / sum(w)
  #   ex2 <- sum((support^2) * w) / sum(w)
  #   return(ex2 - mu^2)
  # } else {
    a_adj <- a - 1L
    log_denom <- .log_denom_truncpois(a_adj, b, lambda)
    mu <- extruncpois(lambda, a, b)
    log_fac_num <- .log_diff(
      stats::ppois(b - 2L, lambda, log.p = TRUE, lower.tail = TRUE),
      stats::ppois(a_adj - 2L, lambda, log.p = TRUE, lower.tail = TRUE)
    )

    ex2 <- exp(2 * log(lambda) + log_fac_num - log_denom) + mu
    return(ex2 - mu^2)
  #}
}
#' Truncated Poisson Median
#'
#' Compute the median of a truncated Poisson distribution
#'
#' Compute the median of a truncated Poisson distribution,
#' \eqn{\text{median}(X \mid a \le X \le b)}, using the formula:
#' \eqn{G^{-1}((G(b, \theta) + G(a-1, \theta))/2, \theta)},
#' where \eqn{G} is the CDF (ppois) and \eqn{G^{-1}} is the quantile function (qpois).
#' With \code{a = 0} and \code{b = Inf} as defaults,
#' the function returns the median of a non-truncated Poisson distribution.
#'
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Upper truncation bound (inclusive). Default is \code{b = Inf}.
#'
#' @return A single numeric value representing the median.
#'
#' @note
#' The median is computed by inverting the average of the CDF values at the truncation boundaries.
#' When \code{lambda < 1}, the median may be close to the lower bound \code{a}.
#'
#' @references
#' Nadarajah, S. and Kotz, S. (2006).
#' \dQuote{R Programs for Computing Truncated Distributions}
#' \emph{Journal of Statistical Software}
#' \url{https://www.jstatsoft.org/v16/c02/}
#'
#' @author Arun Sundar Paulraj and Keefe Murphy
#'
#' @export
#'
#' @examples
#' # Median of a standard truncated Poisson
#' mtruncpois(lambda = 5, a = 0, b = 10)
#'
#' # Median with low intensity: less than 1
#' mtruncpois(lambda = 0.5, a = 0, b = 5)
#'
#' # Median of a zero-truncated Poisson
#' mtruncpois(lambda = 3, a = 1)
mtruncpois <- function(lambda, a = 0L, b = Inf) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)

  # Compute the average of the CDF values at the boundaries
  # G^{-1}((G(b, theta) + G(a-1, theta))/2, theta)
  a_adj <- a - 1L
  p_lower <- ptruncpois(a_adj, lambda, a = a, b = b, lower.tail = TRUE, log.p = FALSE)
  p_upper <- ptruncpois(b, lambda, a = a, b = b, lower.tail = TRUE, log.p = FALSE)
  p_median <- (p_lower + p_upper) / 2

  # Apply the quantile function to get the median
  median <- qtruncpois(p_median, lambda, a = a, b = b, lower.tail = TRUE, log.p = FALSE)
  return(median)
}

#' Truncated Poisson Mode
#'
#' Compute the mode of a truncated Poisson distribution
#'
#' Compute the mode (most probable value) of a truncated Poisson distribution.
#' The untruncated Poisson mode is obtained as \eqn{\text{unique}(c(\lceil\lambda\rceil - 1, \lfloor\lambda\rfloor))},
#' accounting for potential ties. When truncated to \code{[a, b]},
#' the mode is returned as the values clamped to the truncation bounds:
#' \eqn{\text{pmin}(b, \text{pmax}(a, \text{untruncated\_mode}))}, providing all modes
#' that fall within the truncation region.
#'
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Upper truncation bound (inclusive). Default is \code{b = Inf}.
#'
#' @return An integer vector representing all modes in the truncated region.
#'    When multiple modes exist after truncation, all are returned.
#'
#' @note
#' This function uses the closed-form mode formula for the discrete Poisson distribution,
#' accounting for potential ties when \eqn{\lambda} is close to an integer.
#' When \code{lambda < 1}, the untruncated modes are 0, which is also returned
#' (assuming \code{a <= 0}).
#'
#' @references
#' Nadarajah, S. and Kotz, S. (2006).
#' \dQuote{R Programs for Computing Truncated Distributions}
#' \emph{Journal of Statistical Software}
#' \url{https://www.jstatsoft.org/v16/c02/}
#'
#' @author Arun Sundar Paulraj and Keefe Murphy
#'
#' @export
#'
#' @examples
#' # Mode of a standard truncated Poisson
#' mode_truncpois(lambda = 5, a = 0, b = 10)
#'
#' # Mode with potential tie (lambda between two integers)
#' mode_truncpois(lambda = 2.5, a = 0, b = 10)
#'
#' # Mode of a zero-truncated Poisson
#' mode_truncpois(lambda = 3, a = 1)
mode_truncpois <- function(lambda, a = 0L, b = Inf) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)

  # Compute the untruncated mode(s)
  # For Poisson, the mode is unique(c(ceiling(lambda) - 1, floor(lambda)))
  untruncated_modes <- unique(c(ceiling(lambda) - 1, floor(lambda)))

  # Clamp all modes to the truncation bounds [a, b]
  modes <- pmin(b, pmax(a, untruncated_modes))

  # Return unique sorted modes
  return(as.integer(sort(unique(modes))))
}
