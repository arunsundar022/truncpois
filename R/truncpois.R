#' Truncated Poisson Density
#'
#' Compute the probability mass function (PMF) of a truncated Poisson distribution.
#'
#' With \code{a = 0} and \code{b = Inf} as defaults, the function behaves identically to
#' \code{\link[stats]{dpois}}, corresponding to a non-truncated Poisson distribution.
#'
#' @param x Vector of non-negative integers for which probabilities are computed.
#' @param lambda Positive numeric vector, the mean parameter of the Poisson distribution.
#' @param a Non-negative integer lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Strictly positive integer upper truncation bound (inclusive), or \code{Inf}.
#'   Default is \code{b = Inf}.
#' @param log Logical; If \code{log = TRUE}, return log-density. Default is \code{log = FALSE}.
#'
#' @return A numeric vector of densities (or log-densities if \code{log = TRUE})
#'    of the same length as \code{x}. Values of \code{x} outside \code{[a,b]}
#'    return \code{0} or \code{-Inf} if \code{log = TRUE}.
#'
#' @note
#' This implementation adapts the continuous truncated distribution framework of
#' Nadarajah and Kotz (2006) to the discrete Poisson setting.
#' Key differences include:
#' \itemize{
#' \item All internal calculations are performed on the log-scale for numerical stability.
#' \item The normalizing constant is computed as the log-difference of cumulative
#'   probabilities (\code{\link[stats]{ppois}}) evaluated at the truncation bounds,
#'   rather than summing densities.
#' }
#'
#' @seealso
#' \code{\link{ptruncpois}}, \code{\link{qtruncpois}}, \code{\link{rtruncpois}}
#' for the CDF, quantile function, and random generation.
#' \code{\link[stats]{dpois}} for the untruncated Poisson density.
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
#' # Compare standard Poisson and right-truncated Poisson PMF
#' x <- 0:15
#' p_pois <- stats::dpois(x, lambda = 5)
#' p_trunc <- dtruncpois(x, lambda = 5, b = 10)
#'
#' round(rbind(x, p_pois, p_trunc), 3)
#'
#' # Plotting the comparison
#' barplot(rbind(p_pois, p_trunc), beside = TRUE, names.arg = x,
#'         col = c("gray", "blue"),
#'         main = "Standard vs Truncated Poisson PMF",
#'         legend.text = c("Poisson(5)", "TruncPoisson(5, b=10)"))
#'
#' # Log-scale computation for large lambdas and wide bounds to prevent underflow
#' dtruncpois(100, lambda = 100, a = 50, b = 150, log = TRUE)
#'
#' # Check validity: PMF sums to 1 over the entire support
#' upper <- qtruncpois(1 - 1e-12, a = 1, lambda = 6)
#' sum(dtruncpois(seq(1, upper, by = 1), a = 1, lambda = 6))
dtruncpois <- function(x, lambda, a = 0L, b = Inf, log = FALSE) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  if (!is.numeric(x))
    stop("'x' must be numeric", call. = FALSE)
  if (!is.logical(log) || length(log) != 1L)
    stop("'log' must be a single logical", call. = FALSE)

  len <- length(x)
  if (length(lambda) != 1L && length(lambda) != len)
    warning("'lambda' length (", length(lambda), ") not equal to 1 or length(x) (", len, "), recycling may occur", call. = FALSE)
  if (length(a) != 1L && length(a) != len)
    warning("'a' length (", length(a), ") not equal to 1 or length(x) (", len, "), recycling may occur", call. = FALSE)
  if (length(b) != 1L && length(b) != len)
    warning("'b' length (", length(b), ") not equal to 1 or length(x) (", len, "), recycling may occur", call. = FALSE)

  a <- rep(a, length.out = len)
  b <- rep(b, length.out = len)
  lambda <- rep(lambda, length.out = len)

  a_adj <- a - 1L
  log_denom <- .log_denom_truncpois(a_adj, b, lambda)

  dens <- rep(-Inf, len)
  valid <- x > a_adj & x <= b & x == floor(x)
  dens[valid] <- stats::dpois(x[valid], lambda[valid], log = TRUE) - log_denom[valid]

  if (log) return(dens) else return(exp(dens))
}

#' Truncated Poisson Distribution Function
#'
#' Compute the cumulative distribution function (CDF) of a truncated Poisson distribution.
#'
#' With \code{a = 0} and \code{b = Inf} as defaults, the function behaves identically to
#' \code{\link[stats]{ppois}}, corresponding to a non-truncated Poisson distribution.
#'
#' @param q Vector of quantiles.
#' @param lambda Positive numeric vector, the mean parameter of the Poisson distribution.
#' @param a Non-negative integer lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Strictly positive integer upper truncation bound (inclusive), or \code{Inf}.
#'   Default is \code{b = Inf}.
#' @param lower.tail Logical; If \code{TRUE} (default), probabilities are \eqn{P(X \le q)}; otherwise, \eqn{P(X > q)}.
#' @param log.p Logical; If \code{TRUE}, probabilities p are given as log(p). Default is \code{log.p = FALSE}.
#'
#' @return A numeric vector of probabilities (or log-probabilities if \code{log.p = TRUE})
#'    of the same length as \code{q}.
#'
#' @note
#' This implementation adapts the continuous truncated distribution framework of
#' Nadarajah and Kotz (2006) to the discrete Poisson setting.
#' Key differences include:
#' \itemize{
#' \item All internal calculations are performed on the log-scale for numerical stability.
#' \item The normalizing constant is computed as the log-difference of cumulative
#'   probabilities (\code{\link[stats]{ppois}}) evaluated at the truncation bounds.
#' \item All four combinations of \code{lower.tail} and \code{log.p} are correctly
#'   handled entirely on the log-scale, avoiding catastrophic cancellation that
#'   affects naive implementations.
#' }
#'
#' @seealso
#' \code{\link{dtruncpois}}, \code{\link{qtruncpois}}, \code{\link{rtruncpois}}
#' for the PMF, quantile function, and random generation.
#' \code{\link[stats]{ppois}} for the untruncated Poisson CDF.
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
#' # CDF comparison: standard vs truncated Poisson
#' q <- 0:10
#' cdf_pois <- stats::ppois(q, lambda = 4)
#' cdf_trunc <- ptruncpois(q, lambda = 4, b = 7)
#'
#' # The truncated CDF reaches 1.0 exactly at the upper bound b = 7
#' round(rbind(q, cdf_pois, cdf_trunc), 3)
#'
#' # Upper-tail probability P(X > 5 | 2 <= X <= Inf, lambda = 5)
#' ptruncpois(5, lambda = 5, a = 2, lower.tail = FALSE)
#'
#' # Upper-tail log-probability (numerically stable for extreme quantiles)
#' ptruncpois(5, lambda = 5, a = 2, lower.tail = FALSE, log.p = TRUE)
#'
#' # Lower-tail log-probability
#' ptruncpois(10, lambda = 5, a = 2, b = 10, log.p = TRUE)
ptruncpois <- function(q, lambda, a = 0L, b = Inf,
                       lower.tail = TRUE, log.p = FALSE) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  if (!is.numeric(q))
    stop("'q' must be numeric", call. = FALSE)
  if (!is.logical(lower.tail) || length(lower.tail) != 1L)
    stop("'lower.tail' must be a single logical", call. = FALSE)
  if (!is.logical(log.p) || length(log.p) != 1L)
    stop("'log.p' must be a single logical", call. = FALSE)

  len <- length(q)
  if (length(lambda) != 1L && length(lambda) != len)
    warning("'lambda' length (", length(lambda), ") not equal to 1 or length(q) (", len, "), recycling may occur", call. = FALSE)
  if (length(a) != 1L && length(a) != len)
    warning("'a' length (", length(a), ") not equal to 1 or length(q) (", len, "), recycling may occur", call. = FALSE)
  if (length(b) != 1L && length(b) != len)
    warning("'b' length (", length(b), ") not equal to 1 or length(q) (", len, "), recycling may occur", call. = FALSE)

  a <- rep(a, length.out = len)
  b <- rep(b, length.out = len)
  lambda <- rep(lambda, length.out = len)

  a_adj <- a - 1L
  log_denom <- .log_denom_truncpois(a_adj, b, lambda)
  q_floor <- floor(q)

  log_cdf <- rep(-Inf, len)
  above <- q_floor >= b
  mid <- !above & q_floor > a_adj
  log_cdf[above] <- 0
  log_cdf[mid] <- .log_diff(
    stats::ppois(q_floor[mid], lambda[mid], log.p = TRUE, lower.tail = TRUE),
    stats::ppois(a_adj[mid], lambda[mid], log.p = TRUE, lower.tail = TRUE)
  ) - log_denom[mid]

  if (!lower.tail) log_cdf <- .log1mexp(log_cdf)
  if (log.p) return(log_cdf) else return(exp(log_cdf))
}

#' Truncated Poisson Quantile Function
#'
#' Compute the quantile function of a truncated Poisson distribution.
#'
#' With \code{a = 0} and \code{b = Inf} as defaults, the function behaves identically to
#' \code{\link[stats]{qpois}}, corresponding to a non-truncated Poisson distribution.
#'
#' @param p Vector of probabilities.
#' @param lambda Positive numeric vector, the mean parameter of the Poisson distribution.
#' @param a Non-negative integer lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Strictly positive integer upper truncation bound (inclusive), or \code{Inf}.
#'   Default is \code{b = Inf}.
#' @param lower.tail Logical; If \code{TRUE} (default), probabilities are \eqn{P(X \le x)}; otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; If \code{TRUE}, probabilities p are given as log(p). Default is \code{log.p = FALSE}.
#'
#' @return A numeric vector of quantiles of the same length as \code{p}.
#'
#' @note
#' This implementation adapts the continuous truncated distribution framework of
#' Nadarajah and Kotz (2006) to the discrete Poisson setting.
#' Key differences include:
#' \itemize{
#' \item All internal calculations are performed on the log-scale for numerical stability.
#' \item The normalizing constant is computed as the log-difference of cumulative
#'   probabilities (\code{\link[stats]{ppois}}) evaluated at the truncation bounds.
#' \item All four combinations of \code{lower.tail} and \code{log.p} are correctly
#'   handled entirely on the log-scale, consistent with \code{\link[stats]{qpois}}
#'   conventions and avoiding issues present in naive implementations.
#' }
#'
#' @seealso
#' \code{\link{ptruncpois}} for the corresponding CDF.
#' \code{\link{dtruncpois}}, \code{\link{rtruncpois}} for the PMF and random generation.
#' \code{\link[stats]{qpois}} for the untruncated Poisson quantile function.
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
#' # Quantile function for a doubly truncated Poisson distribution
#' p <- seq(0.1, 0.9, by = 0.2)
#' q <- qtruncpois(p, lambda = 10, a = 5, b = 15)
#'
#' # Verify the inverse relationship with ptruncpois:
#' # the CDF at the returned quantiles should be >= the target probabilities
#' cdf_at_q <- ptruncpois(q, lambda = 10, a = 5, b = 15)
#' rbind(Target_p = p, Quantile = q, Actual_CDF = cdf_at_q)
#'
#' # Probabilities supplied on the log-scale (log.p = TRUE)
#' qtruncpois(log(p), lambda = 10, a = 5, b = 15, log.p = TRUE)
#'
#' # Upper-tail quantiles: smallest x with P(X > x) <= p
#' qtruncpois(p, lambda = 10, a = 5, b = 15, lower.tail = FALSE)
#'
#' # Upper-tail quantiles supplied as log-probabilities
#' qtruncpois(log(p), lambda = 10, a = 5, b = 15,
#'            lower.tail = FALSE, log.p = TRUE)
qtruncpois <- function(p, lambda, a = 0L, b = Inf,
                       lower.tail = TRUE, log.p = FALSE) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  if (!is.numeric(p))
    stop("'p' must be numeric", call. = FALSE)
  if (!is.logical(lower.tail) || length(lower.tail) != 1L)
    stop("'lower.tail' must be a single logical", call. = FALSE)
  if (!is.logical(log.p) || length(log.p) != 1L)
    stop("'log.p' must be a single logical", call. = FALSE)

  len <- length(p)
  if (length(lambda) != 1L && length(lambda) != len)
    warning("'lambda' length (", length(lambda), ") not equal to 1 or length(p) (", len, "), recycling may occur", call. = FALSE)
  if (length(a) != 1L && length(a) != len)
    warning("'a' length (", length(a), ") not equal to 1 or length(p) (", len, "), recycling may occur", call. = FALSE)
  if (length(b) != 1L && length(b) != len)
    warning("'b' length (", length(b), ") not equal to 1 or length(p) (", len, "), recycling may occur", call. = FALSE)

  a <- rep(a, length.out = len)
  b <- rep(b, length.out = len)
  lambda <- rep(lambda, length.out = len)

  if (log.p) {
    p <- if (lower.tail) p else .log1mexp(p)
    if (any(p > 0))
      stop("'p' must be <= 0 when 'log.p' is TRUE", call. = FALSE)
  } else {
    if (any(p < 0) || any(p > 1))
      stop("'p' must be in the range [0, 1] when 'log.p' is FALSE",
           call. = FALSE)
    p <- if (lower.tail) log(p) else log1p(-p)
  }

  a_adj <- a - 1L
  pa <- stats::ppois(a_adj, lambda, lower.tail = TRUE, log.p = TRUE)
  pb <- stats::ppois(b, lambda, lower.tail = TRUE, log.p = TRUE)
  p_adj <- .log_sum(pa + .log1mexp(p), p + pb)
  q <- stats::qpois(p_adj, lambda, lower.tail = TRUE, log.p = TRUE)
  return(q)
}

#' Truncated Poisson Random Generation
#'
#' Draw random samples from a truncated Poisson distribution.
#'
#' With \code{a = 0} and \code{b = Inf} as defaults, the function behaves identically to
#' \code{\link[stats]{rpois}}, corresponding to a non-truncated Poisson distribution.
#'
#' @param n Number of random values to return.
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Non-negative integer lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Strictly positive integer upper truncation bound (inclusive), or \code{Inf}.
#'   Default is \code{b = Inf}.
#' @param method Character string specifying the sampling algorithm. One of:
#'   \describe{
#'     \item{\code{"direct"}}{Enumerate the support and sample via
#'       \code{\link[base]{sample}}. Fast when the support is small, but
#'       fails for infinite or very wide supports.}
#'     \item{\code{"inversion"}}{Inversion sampling via \code{\link{qtruncpois}}.
#'       Works for any finite \code{lambda} regardless of support width.}
#'     \item{\code{"bounded"}}{Uniform-CDF inversion through the untruncated
#'       Poisson CDF. Efficient when the truncation region covers most of the
#'       probability mass.}
#'   }
#'   Default is \code{"direct"}.
#'
#' @return A numeric vector of random deviates of length \code{n}.
#'
#' @note
#' Three distinct sampling algorithms are provided to cover different parameter regimes:
#' \itemize{
#' \item \code{"direct"} uses the Gumbel-max trick to sample in one pass from the
#'   enumerated support. It is fast when the support is small but requires finite bounds.
#' \item \code{"inversion"} transforms uniform random variables through
#'   \code{\link{qtruncpois}} and works for any finite \code{lambda} regardless of
#'   support width.
#' \item \code{"bounded"} performs CDF-inversion entirely within the truncation window
#'   via \code{\link[stats]{qpois}}, and is most efficient when truncation removes
#'   little probability mass.
#' }
#' All three methods produce samples confined to \code{[a, b]}.
#'
#' @seealso
#' \code{\link{qtruncpois}} used internally by \code{method = "inversion"}.
#' \code{\link{dtruncpois}}, \code{\link{ptruncpois}} for the PMF and CDF.
#' \code{\link[stats]{rpois}} for the untruncated Poisson sampler.
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
#' # Generate random samples from a zero-truncated Poisson distribution (ZTP)
#' # ZTP is commonly used in count data models where zero is not observable
#' set.seed(123)
#' samples <- rtruncpois(1000, lambda = 2.5, a = 1)
#' table(samples)
#'
#' # Compare sample mean with the theoretical mean
#' mean(samples)
#' extruncpois(lambda = 2.5, a = 1)
#'
#' # Efficiently sample from a heavily truncated distribution using inversion
#' rtruncpois(10, lambda = 100, a = 90, b = 110, method = "inversion")
#'
#' # Use bounded method when truncation bounds enclose most of the probability mass
#' rtruncpois(10, lambda = 50, a = 40, b = 60, method = "bounded")
rtruncpois <- function(n, lambda, a = 0L, b = Inf,
                       method = c("direct", "inversion", "bounded")) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  if (!is.numeric(n) || length(n) != 1L || n < 0L || n != floor(n))
    stop("'n' must be a single non-negative integer", call. = FALSE)
  if (length(lambda) != 1L)
    stop("'lambda' must be a single value", call. = FALSE)
  if (length(a) != 1L)
    stop("'a' must be a single value", call. = FALSE)
  if (length(b) != 1L)
    stop("'b' must be a single value", call. = FALSE)
  if (!is.character(method))
    stop("'method' must be a single character string", call. = FALSE)
  method <- match.arg(method)

  switch(method,

         inversion = {
           x <- qtruncpois(-stats::rexp(n, rate = 1), lambda, a, b,
                           lower.tail = TRUE, log.p = TRUE)
         },

         bounded = {
           lpa <- stats::ppois(a - 1L, lambda, lower.tail = TRUE, log.p = TRUE)
           lpb <- stats::ppois(b,      lambda, lower.tail = TRUE, log.p = TRUE)
           if (is.finite(lpa)) {
             log_u <- lpa + log1p(expm1(lpb - lpa) * stats::runif(n, 0, 1))
           } else {
             log_u <- lpb + log(stats::runif(n))
           }
           x <- stats::qpois(log_u, lambda, lower.tail = TRUE, log.p = TRUE)
         },

         direct = {
           qlo <- stats::qpois(.Machine$double.eps, lambda,
                               lower.tail = TRUE,  log.p = FALSE)
           qhi <- stats::qpois(.Machine$double.eps, lambda,
                               lower.tail = FALSE, log.p = FALSE)
           a_eff   <- max(a, qlo)
           b_eff   <- min(b, qhi)
           if (!is.finite(a_eff) || !is.finite(b_eff))
             stop("Infinite support: try another 'method'", call. = FALSE)
           support <- seq(ceiling(a_eff), floor(b_eff), by = 1L)
           prob    <- dtruncpois(support, lambda, a, b, log = TRUE)
           x       <- support[.gumbel_max(prob, n)]
         }
  )
  return(x)
}
