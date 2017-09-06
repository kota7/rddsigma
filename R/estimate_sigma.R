#' Estimate the Sigma Parameter for Regression Discontinuity Designs with a Mismeasured Running Variable
#'
#' This function provides a high level API for estimating the sigma parameter
#' (the standard deviation of the measurement error in running variables of
#' sharp regression discontinuity designs).
#'
#' @details The method \code{"tsgauss"} estimates the sigma parameter under the
#' assumption that the true running variable and the measurement error are Gaussian.
#' The method \code{"emparam"} relaxes this assumption and allows them to follow
#' some parametric distributions.  Currently, this function supports the Gaussian
#' distribution for the running variable and the Gaussian and Laplace distributions
#' for the measurement error.
#' @param d_vec binary integer vector of treetment assingment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param lower if TRUE, then those with \code{x < cutoff} are assigned (i.e. d = 1).
#' Otherwise those with \code{x > cutoff} are assigned.
#' @param init_sigma initial value of sigma. If NULL (default), randomly chosen
#' @param method character specifying the estimating method.
#' \code{"tsgauss"} and \code{"emparam"} are supported
#' @param x_dist distribution of the true running variable.
#' Used only when \code{method} is \code{"emparam"}.
#' Currently, only "gauss" is supported.
#' @param u_dist distribution of the measurement error.
#' Used only when \code{method} is \code{"emparam"}.
#' Currently, \code{"gauss"} and \code{"lap"} are supported
#' @param ... additional controls. See \code{\link{tsgauss}} or \code{\link{emparam}}
#' @return object of \code{rddsigma} class
#' @export
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- gen_data(500, 0.2, 0)
#'
#' # gaussian-gaussian model
#' estimate_sigma(dat$d, dat$w, 0, method="tsgauss")
#'
#' # em algorithm estimator with parameteric assumptions
#' estimate_sigma(dat$d, dat$w, 0, method="emparam",
#'                x_dist="gauss", u_dist="gauss", verbose=TRUE)
#' estimate_sigma(dat$d, dat$w, 0, method="emparam",
#'                x_dist="gauss", u_dist="lap", verbose=TRUE)
#'
#' # experiment with lower=TRUE
#' dat <- gen_data(500, 0.5, 1, lower=TRUE)
#'
#' estimate_sigma(dat$d, dat$w, 1, lower=TRUE, method="tsgauss")
#'
#' estimate_sigma(dat$d, dat$w, 1, lower=TRUE, method="emparam",
#'                x_dist="gauss", u_dist="gauss", verbose=TRUE)
#' estimate_sigma(dat$d, dat$w, 1, lower=TRUE, method="emparam",
#'                x_dist="gauss", u_dist="lap", verbose=TRUE)
#' }
estimate_sigma <- function(d_vec, w_vec, cutoff, lower=FALSE, init_sigma=NULL,
                           method=c("tsgauss", "emparam"),
                           x_dist=c("gauss"), u_dist=c("gauss", "lap"),
                           ...)
{
  stopifnot(is.character(method))
  stopifnot(length(method) >= 1)

  if (lower) d_vec <- 1-d_vec
  # if lower, d = 1(x < cutoff)
  # so,     1-d = 1(x > cutoff), as the estimating functions assume

  method <- method[1]
  if (method == "tsgauss") {
    out <- tsgauss(d_vec, w_vec, cutoff, init_sigma, ...)
  } else if (method == "emparam") {
    out <- emparam(d_vec, w_vec, cutoff, init_sigma, x_dist, u_dist, ...)
  } else {
    stop(sprintf("method '%s' is not supported", method))
  }

  out
}