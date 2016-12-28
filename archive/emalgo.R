# MLE with full likelihood
# ...maximize all parameters at the same time


# direct likelihood maximization
lfun_i <- function(m_x, s_x, s, D, w, cutoff)
{
  fun <- function(x) dnorm(x, m_x, s_x) * dnorm(w - x, 0, s)

  if (D == 1L) {
    P <- integrate(fun, cutoff, Inf)$value
  } else {
    P <- integrate(fun, -Inf, cutoff)$value
  }

  P <- P*0.99998 + 0.00001  # avoid overflow
  return(log(P))
}


lfun4 <- function(param, dat)
{
  m_x <- param[1]
  s_x <- param[2]
  s   <- param[3]

  out <- 0
  for (i in seq(nrow(dat)))
    out <- out + lfun_i(m_x, s_x, s, dat$D[i], dat$x[i], cutoff)
  return(out)
}

cutoff <- 1
dat <- get_data(200, 0.1, cutoff)
param <- c(1, 1, 1)
optim(param, lfun4, dat = dat, method = "L-BFGS-B",
      lower = c(-10, 0, 0), upper = c(10, 10, 10),
      control = list(fnscale = -1, trace = 1, maxit = 1000))




Q_integrant <- function(x)
{

}

w <- 2
m_x <- 1
s_x <- 1
s <- 1
cutoff <- 1

get_hfuns <- function(w, cutoff, m_x, s_x, s)
{
  h_num <- function(x)
    dnorm(x, m_x, s_x) * dnorm(w - x)
  h1_den <- integrate(h_num, cutoff, Inf)$value
  h2_den <- integrate(h_num, -Inf, cutoff)$value

  h1 <- function(x) h_num(x) / h1_den
  h2 <- function(x) h_num(x) / h2_den
  return(list(h1, h2))
}

hlist <- get_hfuns(w, m_x, s_x, s, cutoff)

hlist[[1]](1)
hlist[[2]](1)

integrate(h, cutoff, Inf)
