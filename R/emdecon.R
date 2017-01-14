#' EM-Deconvolution Estimator
#' @description Estimate the standard deviation of measurement error in the
#' running variable of sharp RDD.  This estimator is constructed under the
#' assumption that measurement error follow a known parametric distribution.
#' @param d_vec binary integer vector of assignment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param u_dist distribution of \eqn{x} and \eqn{u}
#' @param reltol relative tolerance requied
#' @param maxit maximum number of iteration
#' @param verbose if true, progress is reported
#' @param ... currently not used
#' @return object of \code{rddsigma} class
#' @export
#' @examples
#' \dontrun{
#' dat <- gen_data(500, 0.2, 0)
#' emdecon(dat$d, dat$w, 0, u_dist = "gauss", verbose = TRUE)
#' emdecon(dat$d, dat$w, 0, u_dist = "lap", verbose = TRUE)
#' }
emdecon <- function(
  d_vec, w_vec, cutoff, u_dist = c("gauss", "lap"),
  reltol = 1e-5, maxit = 300L, verbose = FALSE, ...)
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

  stopifnot(is.character(u_dist))
  stopifnot(length(u_dist) >= 1)
  u_dist <- u_dist[1]

  ## make wrapper functions for p_x and p_u
  if (u_dist == "gauss") {
    get_px <- function()
    { decon::DeconPdf(w_vec, sigma, error = "normal", fft = TRUE) }
  } else if (u_dist == "lap") {
    get_px <- function()
    { decon::DeconPdf(w_vec, sigma, error = "laplacian", fft = TRUE) }
  } else {
    stop("unsupported u_dist: ", u_dist)
  }

  ## initialize sigma
  sd_w <- sd(w_vec)
  sigma <- sqrt(sd_w^2*0.25)
  convergence <- 1L
  for (it in 1:maxit)
  {
    ## update sigma
    fx <- get_px()
    if (u_dist == "gauss") {
      new_sigma <- emdecon_update_sigma_gauss(
        sigma, d_vec, w_vec, cutoff, fx$x, fx$y)
    } else if (u_dist == "lap") {
      new_sigma <- emdecon_update_sigma_lap(
        sigma, d_vec, w_vec, cutoff, fx$x, fx$y)
    }
    if (verbose) cat("iter", it, ": sigma =", new_sigma, "\n")
    ## check convergence
    if (abs(new_sigma - sigma) < reltol*(abs(sigma) + reltol)) {
      sigma <- new_sigma
      convergence <- 0L
      break
    }
    sigma <- new_sigma
  }

  out <- list(estimate = c(sigma = sigma),
              stderr = NA_real_, avar = matrix(NA_real_, nrow = 1, ncol = 1),
              nobs = n,
              convergence = convergence,
              model = "emdecon", x_dist = "nonparametric", u_dist = u_dist)
  class(out) <- "rddsigma"
  out
}

