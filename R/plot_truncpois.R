#' Plot a Truncated Poisson Distribution
#'
#' Plot the probability mass function (PMF) or cumulative distribution function
#' (CDF) of a truncated Poisson distribution, optionally overlaying the
#' untruncated Poisson for comparison.
#'
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Upper truncation bound (inclusive). Default is \code{b = Inf}.
#' @param type Character string; one of \code{"pmf"} (default) or \code{"cdf"}.
#' @param compare Logical; if \code{TRUE} (default), overlay the untruncated
#'   Poisson distribution for comparison.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{barplot}}
#'   (for \code{type = "pmf"}) or \code{\link[graphics]{plot}} (for \code{type = "cdf"}).
#'
#' @return Invisibly returns a list with elements \code{x} (support values),
#'   \code{truncated} (truncated distribution values), and \code{untruncated}
#'   (untruncated Poisson values, or \code{NULL} if \code{compare = FALSE}).
#'
#' @note
#' The effective support is determined by finding where the untruncated Poisson
#' probability exceeds \code{.Machine$double.eps}, then intersecting with \code{[a, b]}.
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
#' # Plot PMF comparing truncated and untruncated Poisson
#' plot_truncpois(lambda = 5, b = 8)
#'
#' # Plot CDF of a zero-truncated Poisson
#' plot_truncpois(lambda = 3, a = 1, type = "cdf")
#'
#' # Plot PMF without comparison overlay
#' plot_truncpois(lambda = 4, a = 2, b = 9, compare = FALSE)
plot_truncpois <- function(lambda, a = 0L, b = Inf,
                           type = c("pmf", "cdf"),
                           compare = TRUE, ...) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  type <- match.arg(type)

  qlo <- stats::qpois(.Machine$double.eps, lambda, lower.tail = TRUE,  log.p = FALSE)
  qhi <- stats::qpois(.Machine$double.eps, lambda, lower.tail = FALSE, log.p = FALSE)
  x <- seq(max(a, qlo), min(if (is.finite(b)) b else qhi, qhi))

  if (type == "pmf") {
    y_trunc  <- dtruncpois(x, lambda, a = a, b = b)
    y_untrunc <- if (compare) stats::dpois(x, lambda) else NULL

    if (compare) {
      mat <- rbind(Untruncated = y_untrunc, Truncated = y_trunc)
      args <- list(mat, beside = TRUE, names.arg = x,
                   col = c("grey70", "steelblue"),
                   ylab = "Probability", xlab = "x",
                   main = paste0("PMF: TruncPois(lambda=", lambda,
                                 ", a=", a, ", b=", b, ")"),
                   legend.text = c(paste0("Poisson(", lambda, ")"),
                                   paste0("TruncPois(", lambda, ")")))
    } else {
      args <- list(height = y_trunc,
                   col = "steelblue", names.arg = x,
                   ylab = "Probability", xlab = "x",
                   main = paste0("PMF: TruncPois(lambda=", lambda,
                                 ", a=", a, ", b=", b, ")"))
    }
    do.call(graphics::barplot, c(args, list(...)))

  } else {
    y_trunc   <- ptruncpois(x, lambda, a = a, b = b)
    y_untrunc <- if (compare) stats::ppois(x, lambda) else NULL

    plot(x, y_trunc, type = "s", col = "steelblue", lwd = 2,
         ylim = c(0, 1), ylab = "Cumulative Probability", xlab = "x",
         main = paste0("CDF: TruncPois(lambda=", lambda,
                       ", a=", a, ", b=", b, ")"), ...)
    if (compare) {
      graphics::lines(x, y_untrunc, type = "s", col = "grey50",
                      lwd = 2, lty = 2)
      graphics::legend("bottomright",
                       legend = c(paste0("TruncPois(", lambda, ")"),
                                  paste0("Poisson(", lambda, ")")),
                       col = c("steelblue", "grey50"),
                       lwd = 2, lty = c(1, 2), bty = "n")
    }
  }

  invisible(list(x = x, truncated = if (type == "pmf") y_trunc else y_trunc,
                 untruncated = if (compare) y_untrunc else NULL))
}
