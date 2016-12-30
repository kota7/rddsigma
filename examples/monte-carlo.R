library(rddsigma)
library(ggplot2)
library(dplyr)


## todo: parallelize this?

set.seed(123)

results <- NULL
B <- 50
N <- 500
cutoff <- 0
models <- expand.grid(sigma = c(0.3, 1),
                      x_dist = c("gauss", "gaussmix"),
                      u_dist = c("gauss", "laplace"),
                      stringsAsFactors = FALSE)

for (i in 1:nrow(models))
{
  sigma <- models$sigma[i]
  x_dist <- models$x_dist[i]
  u_dist <- models$u_dist[i]
  cat(i, "/", nrow(models), ": sigma =", sigma,
      ", x_dist =", x_dist, ", u_dist =", u_dist, "\n")
  for (b in 1:B)
  {
    cat(sprintf("\r  %4d/%4d", b, B))
    dat <- gen_data(N, sigma, cutoff, u_dist = u_dist, x_dist = x_dist)
    o <- tsgauss(dat$d, dat$w, cutoff)
    tmp <- data.frame(
      data_id = i, sigma = sigma, x_dist = x_dist, method = "two-step gauss",
      param = names(o$estimate), estimate = o$estimate, stderr = o$stderr,
      stringsAsFactors = FALSE)
    results <- rbind(results, tmp)

    o <- em_gauss_lap(dat$d, dat$w, cutoff, quiet = TRUE)
    tmp <- data.frame(
      data_id = i, sigma = sigma, x_dist = x_dist, method = "EM (gauss/lap)",
      param = names(o$estimate), estimate = o$estimate, stderr = o$stderr,
      stringsAsFactors = FALSE)
    results <- rbind(results, tmp)



  }
  cat("\n")
}


# ggplot(results, aes(factor(data_id), estimate, color = method)) +
#   theme_bw() + facet_wrap(~param) +
#   geom_boxplot()

ggplot(subset(results, param == "sigma"),
       aes(factor(data_id), estimate, color = method)) +
  theme_bw() +
  geom_boxplot()

o <- group_by(results, param, data_id, method) %>%
  summarize(sigma = mean(sigma),
            mean_est = mean(estimate),
            sd_est = sd(estimate),
            mean_se = mean(stderr)) %>%
  as.data.frame()

ggplot(subset(o, param == "sigma"), aes(sd_est, mean_se, color = method)) +
  facet_wrap(~ sigma) +
  geom_point()

