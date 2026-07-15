#' truncpois: Truncated Poisson Distribution
#'
#' Provides distribution functions for the left-truncated, right-truncated,
#' and doubly-truncated Poisson distribution: density (\code{\link{dtruncpois}}),
#' cumulative distribution function (\code{\link{ptruncpois}}), quantile function
#' (\code{\link{qtruncpois}}), and random generation (\code{\link{rtruncpois}}).
#' Also includes closed-form functions for the mean (\code{\link{extruncpois}}),
#' variance (\code{\link{vartruncpois}}), median (\code{\link{medtruncpois}}), and
#' mode (\code{\link{modtruncpois}}) of the truncated Poisson distribution, a
#' maximum likelihood estimator (\code{\link{mletruncpois}}), and a plotting
#' function (\code{\link{plottruncpois}}).
#'
#' @details
#' \describe{
#' \item{Type: }{Package}
#' \item{Package: }{\pkg{truncpois}}
#' \item{Version: }{0.1.0}
#' \item{Date: }{2026-04-05}
#' \item{Licence: }{GPL (>= 3)}
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
#' @seealso \url{https://github.com/arunsundar022/truncpois}
#'
#' @aliases truncpois truncpois-package
#' @docType package
#' @keywords package
"_PACKAGE"

.onAttach <- function(lib, pkg) {
  path <- file.path(lib, pkg, "DESCRIPTION")
  version <- read.dcf(path, "Version")
  name <- read.dcf(path, "Package")
  if (interactive()) {
    packageStartupMessage(paste(" _\n| | Truncated Poisson Distributions    _\n| |                                   (_)\n| |_ _ __ _   _ _ __   ___ _ __   ___  _ ___\n| __| '__| | | | '_ \\ / __| '_ \\ / _ \\| / __|\n| |_| |  | |_| | | | |  (_| |_) | (_) | \\__ \\\n\\___|_|   \\__,_|_| |_|\\___| .__/ \\___/|_|___/\n                          | |\n                          |_|   version", version))
  } else {
    packageStartupMessage("\nPackage ", sQuote(name), " version ", version, ".\n")
  }
}
