#' Plot a Truncated Poisson Distribution
#'
#' Plot the probability mass function (PMF), cumulative distribution function
#' (CDF), or quantile function of a truncated Poisson distribution, optionally
#' overlaying the untruncated Poisson for comparison.
#'
#' @param lambda Positive numeric scalar, the mean parameter of the Poisson distribution.
#' @param a Non-negative integer lower truncation bound (inclusive). Default is \code{a = 0}.
#' @param b Strictly positive integer upper truncation bound (inclusive), or \code{Inf}.
#'   Default is \code{b = Inf}.
#' @param type Character string; one of \code{"pmf"} (default), \code{"cdf"}, or \code{"quantile"}.
#' @param compare Logical; if \code{TRUE} (default), overlay the untruncated
#'   Poisson distribution for comparison.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{barplot}}
#'   (for \code{type = "pmf"}) or \code{\link[graphics]{plot}} (for \code{type = "cdf"}
#'   or \code{type = "quantile"}).
#'
#' @return Invisibly returns a list with elements \code{x} (support values, or
#'   probabilities for \code{type = "quantile"}), \code{truncated} (truncated
#'   distribution values), and \code{untruncated} (untruncated Poisson values,
#'   or \code{NULL} if \code{compare = FALSE}).
#'
#' @note
#' When \code{compare = TRUE}, the axis limits are chosen to cover both the
#' truncated and untruncated distributions, so the effect of truncation is
#' always visible. For \code{type = "cdf"}, the x-axis starts at 0 regardless
#' of \code{a}, so that left-truncation is apparent. Internally,
#' \code{plottruncpois()} always calls \code{\link{dtruncpois}},
#' \code{\link{ptruncpois}}, and \code{\link{qtruncpois}} with their
#' \code{log}, \code{log.p}, and \code{lower.tail} arguments left at their
#' defaults (\code{FALSE}, \code{FALSE}, and \code{TRUE} respectively, where
#' applicable); these are not exposed as arguments to this function.
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
#' plottruncpois(lambda = 5, b = 8)
#'
#' # Plot CDF of a zero-truncated Poisson (x-axis starts at 0 to show truncation)
#' plottruncpois(lambda = 3, a = 1, type = "cdf")
#'
#' # Plot the quantile function, comparing truncated and untruncated Poisson
#' plottruncpois(lambda = 5, a = 2, b = 10, type = "quantile")
#'
#' # Plot PMF without comparison overlay
#' plottruncpois(lambda = 4, a = 2, b = 9, compare = FALSE)
plottruncpois <- function(lambda, a = 0L, b = Inf,
                          type = c("pmf", "cdf", "quantile"),
                          compare = TRUE, ...) {
  .check_truncpois_bounds(a, b)
  .check_lambda(lambda)
  type <- match.arg(type)
  dots <- list(...)

  qlo <- stats::qpois(.Machine$double.eps, lambda, lower.tail = TRUE, log.p = FALSE)
  qhi <- stats::qpois(.Machine$double.eps, lambda, lower.tail = FALSE, log.p = FALSE)

  # Support for the truncated distribution
  x_trunc <- seq(max(a, qlo), min(if (is.finite(b)) b else qhi, qhi))

  # When comparing, extend x to cover the untruncated distribution too
  x_full <- seq(qlo, qhi)

  if (type == "pmf") {
    y_trunc <- dtruncpois(x_trunc, lambda, a = a, b = b)
    y_untrunc <- if (compare) stats::dpois(x_full, lambda) else NULL
    default_main <- paste0(
      "PMF: TruncPois(lambda=", lambda,
      ", a=", a, ", b=", b, ")"
    )

    if (compare) {
      y_trunc_full <- dtruncpois(x_full, lambda, a = a, b = b)
      ylim_max <- max(y_trunc_full, y_untrunc, na.rm = TRUE)
      mat <- rbind(Untruncated = y_untrunc, Truncated = y_trunc_full)
      args <- list(mat,
        beside = TRUE, names.arg = x_full,
        col = c("grey70", "steelblue"),
        ylim = c(0, ylim_max * 1.05),
        ylab = "Probability", xlab = "x",
        legend.text = c(
          paste0("Poisson(", lambda, ")"),
          paste0("TruncPois(", lambda, ")")
        )
      )
    } else {
      args <- list(
        height = y_trunc,
        col = "steelblue", names.arg = x_trunc,
        ylab = "Probability", xlab = "x"
      )
    }
    if (is.null(dots$main)) args$main <- default_main
    do.call(graphics::barplot, c(args, dots))

    return(invisible(list(
      x = if (compare) x_full else x_trunc,
      truncated = y_trunc,
      untruncated = if (compare) y_untrunc else NULL
    )))
  } else if (type == "cdf") {
    x_cdf <- seq(0, qhi)
    y_trunc <- ptruncpois(x_cdf, lambda, a = a, b = b)
    y_untrunc <- if (compare) stats::ppois(x_cdf, lambda) else NULL
    default_main <- paste0(
      "CDF: TruncPois(lambda=", lambda,
      ", a=", a, ", b=", b, ")"
    )

    plot_args <- list(x_cdf, y_trunc,
      type = "s", col = "steelblue", lwd = 2,
      xlim = c(0, max(x_cdf)),
      ylim = c(0, 1), ylab = "Cumulative Probability", xlab = "x"
    )
    if (is.null(dots$main)) plot_args$main <- default_main
    do.call(plot, c(plot_args, dots))
    if (compare) {
      graphics::lines(x_cdf, y_untrunc,
        type = "s", col = "grey50",
        lwd = 2, lty = 2
      )
      graphics::legend("bottomright",
        legend = c(
          paste0("TruncPois(", lambda, ")"),
          paste0("Poisson(", lambda, ")")
        ),
        col = c("steelblue", "grey50"),
        lwd = 2, lty = c(1, 2), bty = "n"
      )
    }

    return(invisible(list(
      x = x_cdf,
      truncated = y_trunc,
      untruncated = if (compare) y_untrunc else NULL
    )))
  } else {
    p_grid <- seq(1e-3, 1 - 1e-3, length.out = 200)
    y_trunc <- qtruncpois(p_grid, lambda, a = a, b = b)
    y_untrunc <- if (compare) stats::qpois(p_grid, lambda) else NULL

    default_main <- paste0(
      "Quantile Function: TruncPois(lambda=", lambda,
      ", a=", a, ", b=", b, ")"
    )
    plot_args <- list(p_grid, y_trunc,
      type = "s", col = "steelblue", lwd = 2,
      ylim = range(c(y_trunc, y_untrunc)),
      ylab = "Quantile", xlab = "p"
    )
    if (is.null(dots$main)) plot_args$main <- default_main
    do.call(plot, c(plot_args, dots))
    if (compare) {
      graphics::lines(p_grid, y_untrunc,
        type = "s", col = "grey50",
        lwd = 2, lty = 2
      )
      graphics::legend("bottomright",
        legend = c(
          paste0("TruncPois(", lambda, ")"),
          paste0("Poisson(", lambda, ")")
        ),
        col = c("steelblue", "grey50"),
        lwd = 2, lty = c(1, 2), bty = "n"
      )
    }

    return(invisible(list(
      x = p_grid,
      truncated = y_trunc,
      untruncated = if (compare) y_untrunc else NULL
    )))
  }
}
