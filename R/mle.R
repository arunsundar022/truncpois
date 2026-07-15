#' Maximum Likelihood Estimation for the Truncated Poisson Distribution
#'
#' Estimate the rate parameter \eqn{\lambda} of a truncated Poisson distribution
#' from a sample of observed counts, given known truncation bounds.
#'
#' The log-likelihood \eqn{\sum_i \log f(x_i \mid \lambda, a, b)}, with \eqn{f}
#' the truncated Poisson PMF (\code{\link{dtruncpois}}), is maximized over
#' \eqn{\lambda} using one-dimensional numerical optimization
#' (\code{\link[stats]{optimize}}, which uses Brent's method).
#'
#' @param x Numeric vector of observed integer counts, each within \code{[a, b]}.
#' @param a Non-negative integer lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Strictly positive integer upper truncation bound (inclusive), or \code{Inf}.
#'   Default is \code{b = Inf}.
#' @param interval Numeric vector of length 2 giving the search interval for
#'   \eqn{\lambda}. If \code{NULL} (default), an interval is derived from
#'   \code{mean(x)}.
#'
#' @return A list with elements \code{lambda} (the maximum likelihood estimate),
#'   \code{loglik} (the log-likelihood at that estimate), \code{se} (the standard
#'   error of the estimate, obtained from the observed information via
#'   \code{\link[numDeriv]{hessian}}), and \code{n} (the sample size).
#'
#' @note
#' The estimate is obtained by direct numerical maximization of the exact,
#' closed-form truncated Poisson log-likelihood (via \code{\link{dtruncpois}}),
#' not by method-of-moments or simulation-based approaches. Equivalently,
#' \code{\link[stats]{optimize}} minimizes the negative log-likelihood. The
#' standard error is derived from the curvature of the negative log-likelihood
#' at the estimate: \code{se = sqrt(1 / H)}, where \code{H} is the (scalar)
#' Hessian of the negative log-likelihood evaluated at \code{lambda}, computed
#' via \code{\link[numDeriv]{hessian}}.
#'
#' @seealso
#' \code{\link{dtruncpois}} for the density being maximized.
#' \code{\link{extruncpois}}, \code{\link{medtruncpois}}, \code{\link{modtruncpois}}
#' for the closed-form moments at a given \eqn{\lambda}.
#'
#' @references
#' Nadarajah, S. and Kotz, S. (2006).
#' \dQuote{R Programs for Computing Truncated Distributions}
#' \emph{Journal of Statistical Software}
#' \url{https://www.jstatsoft.org/v16/c02/}
#'
#' @author Arun Sundar Paulraj and Keefe Murphy
#'
#' @importFrom numDeriv hessian
#'
#' @export
#'
#' @examples
#' # Simulate a doubly-truncated Poisson sample and recover lambda
#' lambda_true <- 6
#' a <- 2
#' b <- 15
#' x <- rtruncpois(5000, lambda = lambda_true, a = a, b = b)
#' fit <- mletruncpois(x, a = a, b = b)
#' fit
#'
#' # Compare the fitted lambda to the sample mean-based theoretical mean
#' extruncpois(fit$lambda, a = a, b = b)
#' mean(x)
#'
#' # Compare the empirical proportions to the PMF evaluated at the estimate
#' support   <- a:b
#' tab       <- table(factor(x, levels = support))
#' emp_prop  <- as.numeric(tab) / length(x)
#' theo_prop <- dtruncpois(support, lambda = fit$lambda, a = a, b = b)
#'
#' barplot(rbind(Empirical = emp_prop, Theoretical = theo_prop),
#'   beside = TRUE, names.arg = support,
#'   col = c("grey70", "steelblue"),
#'   ylab = "Probability", xlab = "x",
#'   main = paste0("MLE Recovery: True lambda=", lambda_true,
#'                 ", Fitted lambda=", round(fit$lambda, 3)),
#'   legend.text = c("Empirical", "Theoretical (fitted lambda)"),
#'   args.legend = list(x = "topright", inset = 0.025)
#' )
mletruncpois <- function(x, a = 0L, b = Inf, interval = NULL) {
  .check_truncpois_bounds(a, b)
  if (!is.numeric(x) || length(x) < 1L) {
    stop("'x' must be a non-empty numeric vector", call. = FALSE)
  }
  if (any(x != floor(x))) {
    stop("'x' must contain only integers", call. = FALSE)
  }
  if (any(x < a | x > b)) {
    stop("'x' must lie within [a, b]", call. = FALSE)
  }

  if (is.null(interval)) {
    interval <- c(max(1e-6, mean(x) / 4), mean(x) * 4 + 1)
  }

  negloglik <- function(lambda) -sum(dtruncpois(x, lambda, a = a, b = b, log = TRUE))
  fit <- stats::optimize(negloglik, interval = interval)

  H <- numDeriv::hessian(negloglik, fit$minimum)
  se <- sqrt(1 / H[1, 1])

  list(lambda = fit$minimum, loglik = -fit$objective, se = se, n = length(x))
}
