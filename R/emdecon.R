#' EM-Deconvolution Estimator
#' @export
#' @examples
#' dat <- gen_data(500, 0.2, 0)
#' emdecon(dat$d, dat$w, 0)
emdecon <- function(
  d_vec, w_vec, cutoff,
  u_dist = c("gauss", "lap"),
  reltol = 1e-6, maxit = 1000L,
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

  stopifnot(is.character(u_dist))
  stopifnot(length(u_dist) >= 1)
  u_dist <- u_dist[1]

  ## make wrapper functions for p_x and p_u
  if (u_dist == "gauss") {
    get_px <- function()
    { decon::DeconPdf(w_vec, sigma, error = "normal", fft = TRUE) }
    get_pu <- function()
    { dnorm()}
  } else if (u_dist == "lap") {
    get_px <- function()
    { decon::DeconPdf(w_vec, sigma, error = "laplacian", fft = TRUE) }
  } else {
    stop("unsupported u_dist: ", u_dist)
  }

  ## initialize sigma
  sd_w <- sd(w_vec)
  sigma <- sqrt(sd_w^2*0.25)

  for (it in 1:maxit)
  {
    ## update sigma
    fx <- get_px()

  }
}