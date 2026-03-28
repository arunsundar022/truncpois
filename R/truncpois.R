#' Truncated Poisson density
#'
#' Compute the density of a truncated Poisson distribution.
#'
#' @param x Numeric vector of quantiles.
#' @param lambda Positive numeric scalar, the Poisson mean.
#' @param a Lower truncation bound (inclusive). Default is \code{0}.
#' @param b Upper truncation bound (inclusive). Default is \code{Inf}.
#' @param log Logical; if \code{TRUE}, return the log-density. Default is \code{FALSE}.
#'
#' @return A numeric vector of (log-)densities the same length as \code{x}.
#'   Values outside \code{(a, b]} are assigned density zero (or \code{-Inf} on
#'   the log scale).
#' @export
dtruncpois <- function(x, lambda, a = 0L, b = Inf, log = FALSE) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  if (!is.numeric(x))
    stop("'x' must be numeric", call. = FALSE)
  if (!is.logical(log) || length(log) != 1L)
    stop("'log' must be a single logical", call. = FALSE)

  a           <- a - 1L
  log_denom   <- .log_denom_truncpois(a, b, lambda)

  dens        <- rep(-Inf, length(x))
  valid       <- x > a & x <= b & x == floor(x)
  dens[valid] <- stats::dpois(x[valid], lambda, log = TRUE) - log_denom

  if (log) return(dens) else return(exp(dens))
}

#' Truncated Poisson distribution function
#'
#' Compute the CDF of a truncated Poisson distribution.
#'
#' @param q Numeric vector of quantiles.
#' @param lambda Positive numeric scalar, the Poisson mean.
#' @param a Lower truncation bound (inclusive). Default is \code{0}.
#' @param b Upper truncation bound (inclusive). Default is \code{Inf}.
#' @param lower.tail Logical; if \code{TRUE} (default), return \eqn{P(X \le q)};
#'   otherwise return \eqn{P(X > q)}.
#' @param log.p Logical; if \code{TRUE}, return probabilities on the log scale.
#'   Default is \code{FALSE}.
#'
#' @return A numeric vector of probabilities (or log-probabilities) the same
#'   length as \code{q}.
#' @export
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

  a         <- a - 1L
  log_denom <- .log_denom_truncpois(a, b, lambda)
  q_floor   <- floor(q)

  log_cdf        <- rep(-Inf, length(q))
  above          <- q_floor >= b
  mid            <- !above & q_floor > a
  log_cdf[above] <- 0
  log_cdf[mid]   <- .log_diff(
    stats::ppois(q_floor[mid], lambda, log.p = TRUE, lower.tail = TRUE),
    stats::ppois(a,             lambda, log.p = TRUE, lower.tail = TRUE)
  ) - log_denom

  if (!lower.tail) log_cdf <- .log1mexp(log_cdf)
  if (log.p) return(log_cdf) else return(exp(log_cdf))
}

#' Truncated Poisson quantile function
#'
#' Compute quantiles of a truncated Poisson distribution.
#'
#' @param p Numeric vector of probabilities. When \code{log.p = TRUE}, these
#'   should be log-probabilities (i.e. \eqn{\le 0}).
#' @param lambda Positive numeric scalar, the Poisson mean.
#' @param a Lower truncation bound (inclusive). Default is \code{0}.
#' @param b Upper truncation bound (inclusive). Default is \code{Inf}.
#' @param lower.tail Logical; if \code{TRUE} (default), \code{p} is interpreted
#'   as \eqn{P(X \le x)}; otherwise as \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is interpreted on the log
#'   scale. Default is \code{FALSE}.
#'
#' @return An integer vector of quantiles the same length as \code{p}.
#' @export
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

  a     <- a - 1L
  pa    <- stats::ppois(a, lambda, lower.tail = TRUE, log.p = TRUE)
  pb    <- stats::ppois(b, lambda, lower.tail = TRUE, log.p = TRUE)
  p_adj <- .log_sum(pa + .log1mexp(p), p + pb)
  q     <- stats::qpois(p_adj, lambda, lower.tail = TRUE, log.p = TRUE)
  return(q)
}

#' Truncated Poisson random sampling
#'
#' Draw random samples from a truncated Poisson distribution.
#'
#' @param n Non-negative integer scalar; the number of samples to draw.
#' @param lambda Positive numeric scalar, the Poisson mean.
#' @param a Lower truncation bound (inclusive). Default is \code{0}.
#' @param b Upper truncation bound (inclusive). Default is \code{Inf}.
#' @param method Character string specifying the sampling algorithm. One of:
#'   \describe{
#'     \item{\code{"direct"}}{Enumerate the support and sample via
#'       \code{\link[base]{sample}()}. Fast when the support is small, but
#'       fails for infinite or very wide supports.}
#'     \item{\code{"inversion"}}{Inversion sampling via \code{qtruncpois()}.
#'       Works for any finite \code{lambda} regardless of support width.}
#'     \item{\code{"bounded"}}{Uniform-CDF inversion through the untruncated
#'       Poisson CDF. Efficient when the truncation region covers most of the
#'       probability mass.}
#'   }
#'   Default is \code{"direct"}.
#'
#' @return An integer vector of length \code{n}.
#' @export
rtruncpois <- function(n, lambda, a = 0L, b = Inf,
                       method = c("direct", "inversion", "bounded")) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  if (!is.numeric(n) || length(n) != 1L || n < 0L || n != floor(n))
    stop("'n' must be a single non-negative integer", call. = FALSE)
  if (!is.character(method))
    stop("'method' must be a single character string", call. = FALSE)
  method <- match.arg(method)

  switch(method,

         inversion = {
           x <- qtruncpois(-stats::rexp(n, rate = 1), lambda, a, b,
                           lower.tail = TRUE, log.p = TRUE)
         },

         bounded = {
           pa <- stats::ppois(a - 1L, lambda, lower.tail = TRUE, log.p = FALSE)
           pb <- stats::ppois(b,      lambda, lower.tail = TRUE, log.p = FALSE)
           x  <- stats::qpois(stats::runif(n, pa, pb), lambda, lower.tail = TRUE)
         },

         direct = {
           qlo <- stats::qpois(.Machine$double.eps, lambda,
                               lower.tail = TRUE,  log.p = FALSE)
           qhi <- stats::qpois(.Machine$double.eps, lambda,
                               lower.tail = FALSE, log.p = FALSE)
           a   <- max(a, qlo)
           b   <- min(b, qhi)
           if (!is.finite(a) || !is.finite(b))
             stop("Infinite support: try another 'method'", call. = FALSE)
           support <- seq(ceiling(a), floor(b), by = 1L)
           prob    <- dtruncpois(support, lambda, a, b, log = FALSE)
           x       <- support[sample(length(support), size = n,
                                     prob = prob, replace = TRUE)]
         }
  )
  return(x)
}
