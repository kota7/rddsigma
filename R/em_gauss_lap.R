#' EM Algorithm Estimator for Gaussian-Laplace Model
#' @param d_vec binary integer vector of assignment
#' @param w_vec numeric vector of observed running variable
#' @param cutoff threshold value for assignment
#' @param reltol relative tolerance requied
#' @param maxit maximum number of iteration
#' @param integrate_options controls to pass to \code{\link{integrate}}
#' @param quiet if false, progress is reported
#' @param ... currently not used
#' @export
#' @examples
#' \dontrun{
#' dat <- gen_data(500, 0.2, 0)
#' em_gauss_lap(dat$d, dat$w, 0)
#' }
em_gauss_lap <- function(d_vec, w_vec, cutoff,
                         reltol = 1e-6, maxit = 200L,
                         integrate_options = list(),
                         quiet = FALSE, ...)
{
  ## input validation
  stopifnot(is.numeric(d_vec))
  if (!is.integer(d_vec)) d_vec <- as.integer(d_vec)

  stopifnot(is.numeric(w_vec))

  stopifnot(is.numeric(reltol))
  stopifnot(length(reltol) > 0)
  reltol <- reltol[1]

  stopifnot(is.numeric(maxit))
  stopifnot(length(maxit) > 0)
  maxit <- maxit[1]

  stopifnot(is.logical(quiet))
  stopifnot(length(quiet) > 0)
  quiet <- quiet[1]

  ## integrate options
  ## default values
  integrate_reltol <- 1e-5
  integrate_abstol <- 1e-5
  integrate_subdiv <- 100L
  ## overwrite
  if ("rel.tol" %in% names(integrate_options))
    integrate_reltol <- integrate_options[["rel.tol"]]
  if ("abs.tol" %in% names(integrate_options))
    integrate_abstol <- integrate_options[["abs.tol"]]
  if ("subdivisions" %in% names(integrate_options))
    integrate_subdiv <- integrate_options[["subdivisions"]]


  ## remove NAs, if any
  flg <- !is.na(d_vec) & !is.na(w_vec)
  d_vec <- d_vec[flg]
  w_vec <- w_vec[flg]
  n <- sum(flg)
  stopifnot(n > 1)
  stopifnot(all(d_vec %in% c(0L, 1L)))


  ## mean of x is already identified
  mu_x <- mean(w_vec)

  ## standard error of w is useful later
  sd_w <- sd(w_vec)

  ## set initial guess for parameters
  ## we need: sd_x, sigma
  ## here, we assume equal variance that is consistent with sd_w
  sd_x  <- sqrt(sd_w^2/2)
  sigma <- sd_x

  ## place holders
  hfunc <- NULL
  cur_value <- -Inf
  cur_value_vec <- rep(-Inf, n)

  ## function update value and return the increment
  update_value <- function()
  {
    func <- Map(function(m, w) {
      function(x) dnorm(x, m, sd_x) * bda::dlap(w-x, rate = sqrt(2)/sigma)
    }, mu_x, w_vec)

    values <- rep(list(), n)
    values[d_vec == 1L] <- lapply(func[d_vec == 1L], integrate, cutoff, Inf,
                                  rel.tol = integrate_reltol,
                                  abs.tol = integrate_abstol,
                                  subdivisions = integrate_subdiv)
    values[d_vec == 0L] <- lapply(func[d_vec == 0L], integrate, -Inf, cutoff,
                                  rel.tol = integrate_reltol,
                                  abs.tol = integrate_abstol,
                                  subdivisions = integrate_subdiv)
    values <- unlist(lapply(values, `[[`, "value"))


    new_value <- mean(log(values))
    increment <- new_value - cur_value

    ## update
    cur_value_vec <<- values
    cur_value <<- new_value
    return(increment)
  }

  ## function to update parameters using the current h function
  update_parameters <- function()
  {
    func_sigma <- Map(function(m, w, d) {
      function(x) {
        dnorm(x, m, sd_x) * bda::dlap(w-x, rate = sqrt(2)/sigma) / d * abs(w-x)
      }
    }, mu_x, w_vec, cur_value_vec)
    values_sigma <- rep(list(), n)
    values_sigma[d_vec == 1L] <- lapply(func_sigma[d_vec == 1L], integrate,
                                        cutoff, Inf,
                                        rel.tol = integrate_reltol,
                                        abs.tol = integrate_abstol,
                                        subdivisions = integrate_subdiv)
    values_sigma[d_vec == 0L] <- lapply(func_sigma[d_vec == 0L], integrate,
                                        -Inf, cutoff,
                                        rel.tol = integrate_reltol,
                                        abs.tol = integrate_abstol,
                                        subdivisions = integrate_subdiv)
    new_sigma <- sqrt(2)*mean(unlist(lapply(values_sigma, `[[`, "value")))

    func_sdx <- Map(function(m, w, d) {
      function(x) {
        dnorm(x, m, sd_x) * bda::dlap(w-x, rate = sqrt(2)/sigma) / d * (x-mu_x)^2
      }
    }, mu_x, w_vec, cur_value_vec)
    values_sdx <- rep(list(), n)
    values_sdx[d_vec == 1L] <- lapply(func_sdx[d_vec == 1L], integrate,
                                      cutoff, Inf,
                                      rel.tol = integrate_reltol,
                                      abs.tol = integrate_abstol,
                                      subdivisions = integrate_subdiv)
    values_sdx[d_vec == 0L] <- lapply(func_sdx[d_vec == 0L], integrate,
                                      -Inf, cutoff,
                                      rel.tol = integrate_reltol,
                                      abs.tol = integrate_abstol,
                                      subdivisions = integrate_subdiv)
    new_sdx <- sqrt(mean(unlist(lapply(values_sdx, `[[`, "value"))))

    sigma <<- new_sigma
    sd_x <<- new_sdx
  }

  convergence <- 1L
  for (i in 1:maxit)
  {
    inc <- update_value()
    if (!quiet) {
      cat(sprintf("iter %d: sigma = %.3f, sd_x = %.3f, increment = %1.3e\n",
                  i, sigma, sd_x, inc, "\n"))
    }
    if (abs(inc) < reltol*(abs(cur_value) + reltol)) {
      convergence <- 0L
      break
    }
    update_parameters()
  }

  list(estimate = c(sigma = sigma, mu_x = mu_x, sd_x = sd_x),
       convergence = convergence)
}