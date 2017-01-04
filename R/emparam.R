#' EM Algorithm Estimator for Parametric Model
#' @param d_vec binary integer vector of assignment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param reltol relative tolerance requied
#' @param maxit maximum number of iteration
#' @param integ_method,integ_reltol,integ_depth
#' control parameters for numerical integration
#' @param verbose if true, progress is reported
#' @param ... currently not used
#' @export
#' @examples
#' \dontrun{
#' dat <- gen_data(500, 0.2, 0)
#' emparam(dat$d, dat$w, 0, x_dist = "gauss", u_dist = "lap", verbose = TRUE)
#' }
emparam <- function(
  d_vec, w_vec, cutoff,
  x_dist = c("gauss"), u_dist = c("gauss", "lap"),
  reltol = 1e-6, maxit = 1000L,
  integ_method = c("romberg", "simpson", "trapezoid", "simpson2"),
  integ_reltol = 1e-6, integ_depth = 100L,
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

  if (x_dist == "gauss") {
    if (u_dist == "gauss") {
      cat("to be added!\n")
      return()
    } else if (u_dist == "lap") {
      out <- em_gauss_lap_helper(
        d_vec, w_vec, cutoff,
        reltol, maxit,
        integ_method, integ_reltol, integ_depth, verbose)
    } else {
      stop("unsupported u_dist: ", u_dist)
    }
  } else {
    stop("unsupported x_dist: ", x_dist)
  }

  return(out)
}