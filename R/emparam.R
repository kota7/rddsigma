#' EM Algorithm Estimator for Parametric Model
#' @description Estimate the standard deviation of measurement error in the
#' running variable of sharp RDD.  This estimator is constructed under the
#' assumption that true running variable and measurement error follow
#' known parametric distribution.
#' @param d_vec binary integer vector of assignment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param init_sigma initial values of sigma. If NULL, randomly assigned
#' @param x_dist,u_dist distribution of \eqn{x} and \eqn{u}
#' @param reltol relative tolerance requied
#' @param maxit maximum number of iteration
#' @param integ_method numerical integration method
#' @param integ_reltol relative tolerance for numerical integration
#' @param integ_depth maximum recursion depth for numerical integration
#' @param verbose if true, progress is reported
#' @param ... currently not used
#' @return object of \code{rddsigma} class
#' @export
#' @examples
#' \dontrun{
#' dat <- gen_data(500, 0.2, 0)
#' emparam(dat$d, dat$w, 0, x_dist = "gauss", u_dist = "gauss", verbose = TRUE)
#' emparam(dat$d, dat$w, 0, x_dist = "gauss", u_dist = "lap", verbose = TRUE)
#' }
#' @references
#' Kevin M. Murphy and Robert H. Topel (1985), Estimation and Inference in Two-Step Econometric Models. Journal of Business & Economic Statistics, 3(4), pp.370-379
emparam <- function(
  d_vec, w_vec, cutoff, init_sigma = NULL,
  x_dist = c("gauss"), u_dist = c("gauss", "lap"),
  reltol = 1e-5, maxit = 200L,
  integ_method = c("romberg", "simpson", "trapezoid", "simpson2"),
  integ_reltol = reltol*0.25, integ_depth = 100L,
  verbose = FALSE, ...)
{
  ## input validation
  stopifnot(is.numeric(d_vec))
  if (!is.integer(d_vec)) d_vec <- as.integer(d_vec)
  stopifnot(is.numeric(w_vec))

  ## remove NAs, if any
  flg <- !is.na(d_vec) & !is.na(w_vec)
  d_vec <- d_vec[flg]
  w_vec <- w_vec[flg]
  n <- sum(flg)
  stopifnot(n > 1)
  stopifnot(all(d_vec %in% c(0L, 1L)))

  stopifnot(is.character(x_dist))
  stopifnot(is.character(u_dist))
  stopifnot(length(x_dist) >= 1)
  stopifnot(length(u_dist) >= 1)
  x_dist <- x_dist[1]
  u_dist <- u_dist[1]

  stopifnot(is.numeric(reltol))
  stopifnot(length(reltol) >= 1)
  if (length(reltol) != 1) {
    warning("reltol must be scalar, only the first element is used")
    reltol <- reltol[1]
  }

  stopifnot(is.numeric(maxit))
  stopifnot(length(maxit) >= 1)
  if (length(maxit) != 1) {
    warning("maxit must be scalar, only the first element is used")
    maxit <- maxit[1]
  }
  if (!is.integer(maxit)) maxit <- as.integer(maxit)

  stopifnot(is.character(integ_method))
  stopifnot(length(integ_method) >= 1)
  integ_method <- integ_method[1]
  if (!(integ_method %in% c("romberg", "simpson", "trapezoid", "simpson2")))
    stop("unsupported integ_method: ", integ_method)

  stopifnot(is.numeric(integ_reltol))
  stopifnot(length(integ_reltol) >= 1)
  if (length(integ_reltol) != 1) {
    warning("integ_reltol must be scalar, only the first element is used")
    integ_reltol <- integ_reltol[1]
  }

  stopifnot(is.numeric(integ_depth))
  stopifnot(length(integ_depth) >= 1)
  if (length(integ_depth) != 1) {
    warning("integ_depth must be scalar, only the first element is used")
    integ_depth <- integ_depth[1]
  }
  if (!is.integer(integ_depth)) integ_depth <- as.integer(integ_depth)

  stopifnot(is.logical(verbose))
  stopifnot(length(verbose) >= 1)
  if (length(verbose) != 1) {
    warning("verbose must be scalar, only the first element is used")
    verbose <- verbose[1]
  }


  if (is.null(init_sigma)) {
    init_sigma <- runif(1) * sd(w_vec)
  }


  if (x_dist == "gauss") {
    if (u_dist == "gauss") {
      out <- em_gauss_gauss_helper(
        d_vec, w_vec, cutoff, init_sigma,
        reltol, maxit,
        integ_method, integ_reltol, integ_depth, verbose)
    } else if (u_dist == "lap") {
      out <- em_gauss_lap_helper(
        d_vec, w_vec, cutoff, init_sigma,
        reltol, maxit,
        integ_method, integ_reltol, integ_depth, verbose)
    } else {
      stop("unsupported u_dist: ", u_dist)
    }
  } else {
    stop("unsupported x_dist: ", x_dist)
  }

  out$model <- "emparam"
  out$x_dist <- x_dist
  out$u_dist <- u_dist
  class(out) <- "rddsigma"
  out
}