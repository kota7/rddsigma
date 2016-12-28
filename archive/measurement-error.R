try(setwd("~/Documents/rdd"))
try(setwd("~/../Documents/rdd"))

library(decon)
library(bda)
library(Rcpp)
library(numDeriv)
library(dplyr)
library(tidyr)
library(openxlsx)
library(reshape2)

### data ###
get_data <- function(N, sigma, cutoff = 1,
                     xstar_dist = c("normal", "exp", "normalmix"),
                     error_dist = c("normal", "laplace"), nmix = 10) {
  if (xstar_dist[1] == "normal") {
    x_star <- rnorm(N)
  } else if (xstar_dist[1] == "exp") {
    x_star <- rexp(N)
  } else if (xstar_dist[1] == "normalmix") {
    mu <- seq(-2, 2, length = nmix) %>% sample(N, replace = TRUE)
    x_star <- rnorm(N, mu)
  } else {
    stop("xstar_dist invalid")
  }

  if (error_dist[1] == "laplace") {
    eps <- rlap(N, rate = sqrt(2)/sigma)
  } else if (error_dist[1] == "normal") {
    eps <- rnorm(N, sd = sigma)
  } else {
    stop("eps_dist invalid")
  }
  x <- x_star + eps
  D <- as.integer(x_star > cutoff)

  data.frame(x_star, D, x, eps)
}
###


### MLE (two-step gaussian) ###
lfun <- function(s, dat, mu_x, s_x) {
  x <- dat$x
  D <- dat$D
  mm <- (s/s_x)^2 * (x - mu_x)      # E(eps|x)
  ss <- ((1 - (s/s_x)^2) * s^2)^0.5 # Var(eps|x)^0.5
  P1 <- sapply(seq(N), function(i)
    pnorm(x[i] - cutoff, mean = mm[i], sd = ss,
          log.p = TRUE, lower.tail = TRUE))
  P2 <- sapply(seq(N), function(i)
    pnorm(x[i] - cutoff, mean = mm[i], sd = ss,
          log.p = TRUE, lower.tail = FALSE))
  mean(P1*D + P2*(1-D))
}

lfun_for_hess <- function(param, dat) {
  s    <- param[1]
  mu_x <- param[2]
  s_x  <- param[3]
  lfun(s, dat, mu_x, s_x)
}

lfun_for_each <- function(s, dat, mu_x, s_x) {
  x <- dat$x
  D <- dat$D
  mm <- (s/s_x)^2 * (x - mu_x)      # E(eps|x)
  ss <- ((1 - (s/s_x)^2) * s^2)^0.5 # Var(eps|x)^0.5
  P1 <- sapply(seq(N), function(i)
    pnorm(x[i] - cutoff, mean = mm[i], sd = ss,
          log.p = TRUE, lower.tail = TRUE))
  P2 <- sapply(seq(N), function(i)
    pnorm(x[i] - cutoff, mean = mm[i], sd = ss,
          log.p = TRUE, lower.tail = FALSE))
  P1*D + P2*(1-D)
}

get_avar <- function(s, dat, mu_x, s_x) {
  R1 <- -diag(1/c(s_x^2, s_x^2/2))
  tmp <- numDeriv::hessian(lfun_for_hess, c(s, mu_x, s_x), dat = dat)
  R2 <- -tmp[1,1]
  R3 <- -tmp[1, c(2,3)]
  tmp1 <- cbind((dat$x-mu_x)/s_x^2, -1/s_x + (dat$x-mu_x)^2/s_x^3)
  tmp2 <- numDeriv::jacobian(lfun_for_each, s, dat = dat, mu_x = mu_x, s_x = s_x)
  R4 <- crossprod(tmp1, tmp2) / nrow(dat)

  solve(R2) +
    solve(R2) %*% (t(R3) %*% solve(R1) %*% R3 -
                   t(R4) %*% solve(R1) %*% R3 -
                   t(R3) %*% solve(R1) %*% R4) %*% solve(R2)
}
###



### deconvolution estimator ###

computeP_normal <- function(x, dx, w, dw, s, cutoff) {
  N <- length(w)
  P <- numeric(N)
  for (i in seq(N)) {
    dxw <- dnorm(w[i] - x, sd = s) * dx / dw[i]
    P[i] <- sum(dxw[x > cutoff]) / sum(dxw)
  }
  P*0.99998 + 0.00001  # avoid overflow
}

computeP_laplace <- function(x, dx, w, dw, s, cutoff) {
  N <- length(w)
  P <- numeric(N)
  for (i in seq(N)) {
    dxw <- dlap(w[i] - x, rate = sqrt(2)/s) * dx / dw[i]
    P[i] <- sum(dxw[x > cutoff]) / sum(dxw)
  }
  P*0.99998 + 0.00001  # avoid overflow
}

cppFunction("
  NumericVector computeP_normalC(NumericVector x, NumericVector dx,
                                 NumericVector w, NumericVector dw,
                                 double s, double cutoff) {
    int N = w.size();
    NumericVector P(N);
    NumericVector dwx(x.size());
    NumericVector dwx_ss(sum(x>0));

    for (int i = 0; i < N; i++) {
      dwx = dnorm(w[i] - x, 0, s) * dx / dw[i];
      dwx_ss = dwx[x > cutoff];
      P[i] = sum(dwx_ss) / sum(dwx);
    }
    return P * 0.99998 + 0.00001;
  }
")

cppFunction("
  NumericVector computeP_laplaceC(NumericVector x, NumericVector dx,
                                  NumericVector w, NumericVector dw,
                                  double s, double cutoff) {
    int N = w.size();
    NumericVector P(N);
    NumericVector dwx(x.size());
    NumericVector dwx_ss(sum(x>0));

    for (int i = 0; i < N; i++) {
      dwx = exp(-abs(w[i] - x)/(s/sqrt(2)))/(s*sqrt(2)) * dx / dw[i];
      dwx_ss = dwx[x > cutoff];
      P[i] = sum(dwx_ss) / sum(dwx);
    }
    return P * 0.99998 + 0.00001;
  }
")

# check functions
s <- 0.5
cutoff <- 1
dat <- get_data(N = 1000, sigma = s, cutoff = cutoff)
fx <- DeconPdf(dat$x, s, error = "laplacian", fft = TRUE)
fw <- approx(density(dat$x), NULL, dat$x)

all.equal(computeP_normal(fx$x, fx$y, fw$x, fw$y, s, cutoff),
          computeP_normalC(fx$x, fx$y, fw$x, fw$y, s, cutoff))
all.equal(computeP_laplace(fx$x, fx$y, fw$x, fw$y, s, cutoff),
          computeP_laplaceC(fx$x, fx$y, fw$x, fw$y, s, cutoff))


# likelihood functions for deconvolution estimator
lfun2 <- function(s, dat) {
  # normal error
  w <- dat$x
  D <- dat$D
  N <- length(w)

  fx <- DeconPdf(w, s, error = "normal", fft = TRUE)
  fw <- approx(density(w), NULL, w)
  #P <- computeP_normal(fx$x, fx$y, fw$x, fw$y, s)
  P <- computeP_normalC(fx$x, fx$y, fw$x, fw$y, s, cutoff)
  mean(log(P)*D + log(1-P)*(1-D))
}


lfun3 <- function(s, dat) {
  # laplace error
  w <- dat$x
  D <- dat$D
  N <- length(w)

  fx <- DeconPdf(w, s, error = "laplacian", fft = TRUE)
  fw <- approx(density(w), NULL, w)
  #P <- computeP_laplace(fx$x, fx$y, fw$x, fw$y, s)
  P <- computeP_laplaceC(fx$x, fx$y, fw$x, fw$y, s, cutoff)
  mean(log(P)*D + log(1-P)*(1-D))
}
###


### EM algorithm ###

###




### testing ###
N <- 1000
sigma <- 0.2
cutoff <- 1
dat <- get_data(N, sigma, cutoff,
                xstar_dist = "normal", error_dist = "normal")
with(dat, table(D, as.integer(x > cutoff)))
with(dat, mean(D == as.integer(x > cutoff)))


cat("* gaussian / gaussian *\n")
cat("true param", sigma, "\n")
mu_x <- with(dat, mean(x))
s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
avar <- get_avar(o$par, dat = dat, mu_x = mu_x, s_x = s_x)
se <- (avar/N)^0.5
cat("two-step gaussian", o$par, "(", se, ")", "\n")


# # bootstrap
# B <- 300
# tmp <- data.frame(s = numeric(B), cil = numeric(B), ciu = numeric(B))
# for (i in seq(nrow(tmp))) {
#   cat(sprintf("\r%4d", i))
#   dat <- get_data(N, sigma, cutoff,
#                   xstar_dist = "normal", error_dist = "normal")
#   mu_x <- with(dat, mean(x))
#   s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
#   o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
#              lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
#              control = list(fnscale = - 1, reltol = 1e-4))
#   avar <- get_avar(o$par, dat = dat, mu_x = mu_x, s_x = s_x)
#   se <- (avar/N)^0.5
#   tmp$s[i] <- o$par
#   tmp$cil[i] <- o$par - 1.96*se
#   tmp$ciu[i] <- o$par + 1.96*se
# }
# mean(0.2 > tmp$cil & 0.2 < tmp$ciu)


# Standard error <-- check two-step estimator (Murphy and Topel, 1985)

o <- optim(s_x/2, lfun2, dat = dat,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, trace = 1, reltol = 1e-4))
cat("deconvolution", o$par, "\n")




N <- 1000
sigma <- 0.2
cutoff <- 1
dat <- get_data(N, sigma, cutoff,
                xstar_dist = "normal", error_dist = "laplace")
with(dat, table(D, as.integer(x > cutoff)))
with(dat, mean(D == as.integer(x > cutoff)))

mu_x <- with(dat, mean(x))
s_x  <- with(dat, mean((x - mu_x)^2)^0.5)

cat("* gaussian / laplace *\n")
cat("true param", sigma, "\n")
o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
           lower = 0, upper = s_x, method = "Brent",
           hessian = TRUE, control = list(fnscale = - 1, reltol = 1e-4))
avar <- get_avar(o$par, dat = dat, mu_x = mu_x, s_x = s_x)
se <- (avar/N)^0.5
cat("two-step gaussian", o$par, "(", se, ")", "\n")


o <- optim(s_x/2, lfun3, dat = dat,
           lower = 0, upper = s_x, method = "Brent",
           hessian = TRUE, control = list(fnscale = - 1, reltol = 1e-4))

cat("deconvolution", o$par, "\n")




N <- 1000
sigma <- 0.2
cutoff <- 1
dat <- get_data(N, sigma, cutoff,
                xstar_dist = "exp", error_dist = "normal")
with(dat, table(D, as.integer(x > cutoff)))
with(dat, mean(D == as.integer(x > cutoff)))


cat("* exponential / gaussian *\n")
cat("true param", sigma, "\n")
mu_x <- with(dat, mean(x))
s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
se <- (avar/N)^0.5
cat("two-step gaussian", o$par, "(", se, ")", "\n")

o <- optim(s_x/2, lfun3, dat = dat,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
cat("deconvolution", o$par, "\n")





N <- 1000
sigma <- 0.2
cutoff <- 1
dat <- get_data(N, sigma, cutoff,
                xstar_dist = "exp", error_dist = "laplace")
with(dat, table(D, as.integer(x > cutoff)))
with(dat, mean(D == as.integer(x > cutoff)))


cat("* exponential / laplace *\n")
cat("true param", sigma, "\n")
mu_x <- with(dat, mean(x))
s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
se <- (avar/N)^0.5
cat("two-step gaussian", o$par, "(", se, ")", "\n")

o <- optim(s_x/2, lfun3, dat = dat,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
cat("deconvolution", o$par, "\n")





N <- 1000
sigma <- 0.2
cutoff <- 1
dat <- get_data(N, sigma, cutoff,
                xstar_dist = "normalmix", error_dist = "normal")
with(dat, table(D, as.integer(x > cutoff)))
with(dat, mean(D == as.integer(x > cutoff)))


cat("* gaussian mix / gaussian *\n")
cat("true param", sigma, "\n")
mu_x <- with(dat, mean(x))
s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
se <- (avar/N)^0.5
cat("two-step gaussian", o$par, "(", se, ")", "\n")

o <- optim(s_x/2, lfun3, dat = dat,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
cat("deconvolution", o$par, "\n")



N <- 1000
sigma <- 0.2
cutoff <- 1
dat <- get_data(N, sigma, cutoff,
                xstar_dist = "normalmix", error_dist = "laplace")
with(dat, table(D, as.integer(x > cutoff)))
with(dat, mean(D == as.integer(x > cutoff)))


cat("* gaussian mix / laplace *\n")
cat("true param", sigma, "\n")
mu_x <- with(dat, mean(x))
s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
o <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
se <- (avar/N)^0.5
cat("two-step gaussian", o$par, "(", se, ")", "\n")

o <- optim(s_x/2, lfun3, dat = dat,
           lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
           control = list(fnscale = - 1, reltol = 1e-4))
cat("deconvolution", o$par, "\n")
###


stop()
### Monte Carlo exercise ###
cutoff <- 1
N <- 1000
B <- 200
seed <- 87
set.seed(seed)
out <- NULL
for (sigma in c(0.2, 0.8)) {
  for (xstar_dist in c("normalmix", "normal", "exp")) {
    for (error_dist in c("normal", "laplace")) {
      cat("\n*", sigma, xstar_dist, error_dist, "\n")
      for (b in seq(B)) {
        cat(sprintf("\r%4d/%4d  ", b, B))
        dat <- get_data(N, sigma, cutoff, xstar_dist, error_dist)
        # two step gaussian
        mu_x <- with(dat, mean(x))
        s_x  <- with(dat, mean((x - mu_x)^2)^0.5)
        o1 <- optim(s_x/2, lfun, dat = dat, mu_x = mu_x, s_x = s_x,
                    lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
                    control = list(fnscale = - 1, reltol = 1e-4))

        # deconvolution with normal
        o2 <- optim(s_x/2, lfun2, dat = dat,
                    lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
                    control = list(fnscale = - 1, trace = 1, reltol = 1e-4))

        # deconvolution with laplace
        o3 <- optim(s_x/2, lfun3, dat = dat,
                    lower = 0, upper = s_x, method = "Brent", hessian = TRUE,
                    control = list(fnscale = - 1, trace = 1, reltol = 1e-4))

        tmp <- data.frame(b, sigma, xstar_dist, error_dist,
          ts_gauss = o1$par, decon_gauss = o2$par, decon_lap = o3$par)
        out <- if (is.null(out)) tmp else rbind(out, tmp)
      }
    }
  }
}

write.csv(out, sprintf("output/montecarlo%d-%d-%d.csv", B, seed, N), row.names = FALSE)

tmp1 <- out %>% melt(id.var = c("b", "sigma", "xstar_dist", "error_dist"),
                     variable.name = "est_method", value.name = "estimate") %>%
  group_by(sigma, xstar_dist, error_dist, est_method) %>%
  summarize(mean(estimate), sd(estimate)) %>%
  mutate(cil = `mean(estimate)` - 1.96 * `sd(estimate)`,
         ciu = `mean(estimate)` + 1.96 * `sd(estimate)`) %>%
  mutate(true_dist = paste(xstar_dist, error_dist)) %>% ungroup %>%
  mutate(est_str = sprintf("%.3f", `mean(estimate)`),
         sd_str = sprintf("(%.3f)", `sd(estimate)`),
         ci_str = sprintf("[%.3f, %.3f]", cil, ciu))

tmp2 <- tmp1 %>%
  select(-xstar_dist, -error_dist) %>%
  # select(sigma, true_dist, est_method, ends_with("_str")) %>%
  melt(id.var = c("sigma", "true_dist", "est_method"),
       variable.name = "stats", value.name = "value")

tmp3 <- tmp2 %>% filter(regexpr("_str$", stats) > 0) %>%
  spread(key = "true_dist", value = "value")
tmp4 <- tmp2 %>% filter(regexpr("_str$", stats) > 0) %>%
  spread(key = "est_method", value = "value")

xlsout <- list(raw = out, long = tmp2,
               wide_stats = tmp1, wide_dist = tmp3, wide_method = tmp4)
openxlsx::write.xlsx(
  xlsout, sprintf("output/montecarlo%d-%d-%d.xlsx", B, seed, N))
###



# benchmark
# library(microbenchmark)
# microbenchmark(
#   computeP_normal(fx$x, fx$y, fw$x, fw$y, s),
#   computeP_normalC(fx$x, fx$y, fw$x, fw$y, s),
#   computeP_laplace(fx$x, fx$y, fw$x, fw$y, s),
#   computeP_laplaceC(fx$x, fx$y, fw$x, fw$y, s)
# )

#
# ### STAN ###
# # not working very well
# library(rstan)
#
# standata <- with(dat, list(N = length(x), x = x, D = D))
#
# pars <- c("mu_x", "sd_x", "sd_e", "sd_sum")
# warmup <- 5000
# iter   <- warmup + 10000
# seed   <- 87
# digit  <- 3
# chaing <- 1
# fit <- stan(file = "stan/gaussian.stan", model_name = "gaussian",
#             data = standata, pars = pars, seed = seed, chains = chaing,
#             warmup = warmup, iter = iter, verbose = FALSE)
#
# print(fit, pars = pars, digits_summary = digit)
# plot(fit, pars = pars)
# traceplot(fit, inc_warmup=TRUE, pars=pars)


