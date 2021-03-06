#' Two Step Gaussian Estimator
#' @description Estimate the standard deviation of measurement error in the
#' running variable of sharp RDD.  This estimator is constructed under
#' the assumption that true running variable and measurement error
#' follow the Gaussian distribution.
#' @param d_vec binary integer vector of assignment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param init_sigma initial value of sigma. If NULL, randomly assigned
#' @param ... additional controls for \code{optim}
#' @return object of \code{rddsigma} class
#' @export
#' @examples
#' set.seed(123)
#' dat <- gen_data(500, 0.2, 0)
#' tsgauss(dat$d, dat$w, 0)
#' @references
#' Kevin M. Murphy and Robert H. Topel (1985), Estimation and Inference in Two-Step Econometric Models. Journal of Business & Economic Statistics, 3(4), pp.370-379
tsgauss <- function(d_vec, w_vec, cutoff, init_sigma = NULL, ...)
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

  ## ML estimate for mu_x (= mu_w) and sigma_w
  mu_x = mean(w_vec)
  sd_w = sd(w_vec) * (n-1) / n  # ML estimate should divides by n

  lfunc <- function(sigma)
  {
    # returns: mean log-likelihood

    ## distribution of f(U | W)
    mu_u <- sigma^2/sd_w^2 * (w_vec - mu_x)
    sd_u <- sqrt((1 - sigma^2/sd_w^2) * sigma^2)
    ## use helper function
    tsgauss_lfunc_helper(d_vec, w_vec, cutoff, mu_u, sd_u, sigma)
  }


  lfunc_para3 <- function(param)
  {
    # this function takes sigma, mu_x, sd_w as argument
    # used for hessian evaluation
    sigma <- param[1]
    mu_x <- param[2]
    sd_w <- param[3]

    mu_u <- sigma^2/sd_w^2 * (w_vec - mu_x)
    sd_u <- sqrt((1 - sigma^2/sd_w^2) * sigma^2)
    tsgauss_lfunc_helper(d_vec, w_vec, cutoff, mu_u, sd_u, sigma)
  }
  lfunc_each <- function(param)
  {
    # returns the full vector of likelihood
    # used for jacobian computation
    sigma <- param[1]
    mu_x <- param[2]
    sd_w <- param[3]

    ## distribution of f(U | W)
    mu_u <- sigma^2/sd_w^2 * (w_vec - mu_x)
    sd_u <- sqrt((1 - sigma^2/sd_w^2) * sigma^2)
    ## use helper function
    tsgauss_lfunc_each_helper(d_vec, w_vec, cutoff, mu_u, sd_u, sigma)
  }

  get_avar <- function(sigma)
  {
    # implements Murphy and Topel (1985), section 5.1
    R1 <- diag(1/c(sd_w^2, sd_w^2/2))
    #tmp <- numDeriv::hessian(lfunc_para3, c(sigma, mu_x, sd_w))
    #R2 <- -tmp[1,1]
    #R3 <- -tmp[c(2,3), 1, drop = FALSE]
    #tmp2 <- numDeriv::jacobian(lfunc_each, sigma)
    #R4 <- crossprod(tmp1, tmp2) / n
    jaco1 <- cbind((w_vec-mu_x)/sd_w^2, -1/sd_w + (w_vec-mu_x)^2/sd_w^3)
    jaco2 <- numDeriv::jacobian(lfunc_each, c(sigma, mu_x, sd_w))
    R2 <- crossprod(jaco2[, 1, drop = FALSE]) / n
    R3 <- crossprod(jaco2[, c(2,3)], jaco2[, 1, drop = FALSE]) / n
    R4 <- crossprod(jaco1, jaco2[, 1, drop = FALSE]) / n

    o11 <- solve(R1)
    o12 <- solve(R1) %*% (R4 - R3) %*% solve(R2)
    o22 <- solve(R2) + solve(R2) %*% (t(R3) %*% solve(R1) %*% R3 -
                                        t(R4) %*% solve(R1) %*% R3 -
                                        t(R3) %*% solve(R1) %*% R4) %*% solve(R2)
    o <- rbind(cbind(o11, o12), cbind(t(o12), o22))

    # compute the avar for sd_x by delta method
    sd_x <- sqrt(sd_w^2 - sigma^2)  # point estimate

    # jacobian matrix for conversion:
    # (mu_x, sd_w, sigma) -> (mu_x, sd_w, sigma, sd_x)
    Gmat <- rbind(diag(3),
                  matrix(c(0, sd_w/sd_x, -sigma/sd_x), nrow = 1, ncol = 3))
    o_aug <- Gmat %*% o %*% t(Gmat)
    rownames(o_aug) <- c("mu_x", "sd_w", "sigma", "sd_x")
    colnames(o_aug) <- c("mu_x", "sd_w", "sigma", "sd_x")


    # reorder, put sigma in front
    o_aug[c("sigma", "mu_x", "sd_x", "sd_w"),
          c("sigma", "mu_x", "sd_x", "sd_w")]
  }


  ## range is between 1e-8 to sd_w
  ## do not set to zero to avoid numerical error
  if (is.null(init_sigma)) {
    init_sigma <- runif(1) * sd(w_vec)
  }
  o <- optim(init_sigma, lfunc, lower = 1e-8, upper = sd_w,
             method = "Brent", hessian = TRUE,
             control = list(fnscale = -1, ...))
  avar <- get_avar(o$par)

  out <- list(estimate = c(sigma = o$par, mu_x = mu_x,
                           sd_x = sqrt(sd_w^2 - o$par^2), sd_w = sd_w),
              stderr = sqrt(diag(avar)/n), avar = avar, nobs = n,
              value = o$value,
              convergence = o$convergence,
              model = "tsgauss", x_dist = "gauss", u_dist = "gauss")
  class(out) <- "rddsigma"
  out
}



## pure R implementation, significantly slower than that using c++ helper
## kept for reference and debugging
tsgauss_r <- function(d_vec, w_vec, cutoff, init_ratio = 0.25, ...)
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

  ## range is between 1e-8 to sd_w
  ## do not set to zero to avoid numerical error
  initial_sigma <- sqrt(sd_w^2*init_ratio)
  o <- optim(initial_sigma, lfunc, lower = 1e-8, upper = sd_w,
             method = "Brent", hessian = TRUE,
             control = list(fnscale = -1, ...))
  return(o)
}

