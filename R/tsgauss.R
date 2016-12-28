#' Two Step Gaussian Model for Measurement Error Variance for a Sharp RDD
#' @param d_vec binary integer vector of assignment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param ... additional controls for \code{optim}
#' @return estimate for the standard deviation
#' @export
#' @examples
#' dat <- gen_data(1000, 0.2, 0)
#' tsgauss(dat$d, dat$w, 0)
tsgauss <- function(d_vec, w_vec, cutoff, ...)
{
  ## remove NAs, if any
  flg <- !is.na(d_vec) & !is.na(w_vec)
  d_vec <- d_vec[flg]
  w_vec <- w_vec[flg]
  n <- sum(flg)
  stopifnot(n > 1)

  ## ML estimate for mu_x (= mu_w) and sigma_w
  mu_x = mean(w_vec)
  sd_w = sd(w_vec) * (n-1) / n  # ML estimate should divides by n

  lfunc <- function(sigma)
  {
    ## distribution of f(U | W)
    mu_u <- sigma^2/sd_w^2 * (w_vec - mu_x)
    sd_u <- sqrt((1 - sigma^2/sd_w^2) * sigma^2)
    P1 <- Map(pnorm, w_vec - cutoff,
              mean = mu_u, sd = sd_u, lower.tail = TRUE, log.p = TRUE)
    P2 <- Map(pnorm, w_vec - cutoff,
              mean = mu_u, sd = sd_u, lower.tail = FALSE, log.p = TRUE)
    mean(unlist(P1)*d_vec + unlist(P2)*(1-d_vec))
  }

  ## initial value is set to sd_w/2
  ## range is between 0 to sd_w
  o <- optim(sd_w/2, lfunc, lower = 0, upper = sd_w,
             method = "Brent", hessian = TRUE,
             control = list(fnscale = -1, ...))
  return(o)



  ## todo: compute standard error

}

