#' Generate Simulation Data for a Sharp RDD with Mismeasured Running Variable
#' @param n sample size
#' @param sigma standard deviation of the measurement error term
#' @param cutoff threshold value for the assignment
#' @param u_dist distribution of the measurement error.
#' Only "gauss" or "lap" is supported.
#' @param x_dist distribution of the true running variable.
#' Either "gauss", "gaussmix" or any function that a generate random sample.
#' @param mu_x,sd_x mean and standard deviation of gaussian distribution.
#' Used if \code{x_dist} is 'gauss'
#' @param rate_x rate parameter used when \code{x_dist} is 'exp'
#' @param min_x,max_x min and max parameter used when \code{x_dist} is 'unif'
#' @param nmix,mu_min,mu_max,sds  parameters used
#' when \code{x_dist} is 'gaussmix'.
#' \code{nmix} is the number of mixture. \code{mu_min} and \code{mu_max} is
#' the range of means. \code{sds} is the standard deviations.
#' @param ... additional argument to be passed to the generator function for
#' running variable
#' @return list of 4 numeric vectors of length n.
#' \itemize{
#' \item{d}: assignment variable
#' \item{w}: observed running variable
#' \item{x}: true running variable
#' \item{u}: measurement error
#' }
#' @export
#' @examples
#' gen_data(100, 0.2, 0)
#' gen_data(100, 0.2, 1, u_dist = "lap", x_dist = "exp")
gen_data <- function(n, sigma, cutoff,
                     u_dist = c("gauss", "lap"),
                     x_dist = c("gauss", "exp", "unif", "gaussmix"),
                     mu_x = 0, sd_x = 1,
                     rate_x = 1,
                     min_x = -sqrt(3), max_x = sqrt(3),
                     nmix = 5, mu_min = -1.5, mu_max = 1.5, sds = 0.25,
                     ...)
{
  ## validation for numeric arguments
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(sigma))
  stopifnot(is.numeric(cutoff))

  stopifnot(length(n) > 0)
  stopifnot(length(sigma) > 0)
  stopifnot(length(cutoff) > 0)

  if (length(n) > 1) {
    warning("n must be a scalar. only the first element is used")
    n <- n[1]
  }
  if (length(sigma) > 1) {
    warning("sigma must be a scalar. only the first element is used")
    sigma <- sigma[1]
  }
  if (length(cutoff) > 1) {
    warning("cutoff must be a scalar. only the first element is used")
    cutoff <- cutoff[1]
  }
  stopifnot(n > 0)
  stopifnot(sigma > 0)


  ## validation for generation function specification
  stopifnot(is.character(u_dist))
  stopifnot(length(u_dist) > 0)
  u_dist <- u_dist[1]
  if (u_dist == "gauss") {
    u_generator <- function(n) rnorm(n, sd = sigma)
  } else if (u_dist == "lap") {
    u_generator <- function(n) bda::rlap(n, rate = sqrt(2)/sigma)
  } else {
    stop("u_dist must be either 'gauss' or 'laplace'")
  }

  if (is.character(x_dist)) {
    stopifnot(length(x_dist) > 0)
    x_dist <- x_dist[1]
    if (x_dist == "gauss") {
      x_generator <- function(n) rnorm(n, mean = mu_x, sd = sd_x)
    } else if (x_dist == "gaussmix") {
      x_generator <- function(n) {
        mixtools::rnormmix(n,
                           lambda = rep(1/nmix, nmix),
                           mu = seq(mu_min, mu_max, length = nmix),
                           sigma = rep_len(sds, nmix))
      }
    } else if (x_dist == "exp") {
      x_generator <- function(n) rexp(n, rate = rate_x)
    } else if (x_dist == "unif") {
      x_generator <- function(n) runif(n, min = min_x, max = max_x)
    } else {
      stop("unknown character string for x_dist", x_dist)
    }
  } else if (is.function(x_dist)) {
    x_generator <- function(n) x_dist(n, ...)
  } else {
    stop("invalid input for x_dist")
  }

  u <- u_generator(n)
  x <- x_generator(n)
  d <- as.integer(x > cutoff)
  w <- x + u

  list(d = d, w = w, x = x, u = u)
}

