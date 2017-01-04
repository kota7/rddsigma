library(rddsigma)
library(ggplot2)
library(dplyr)



set.seed(123)

results <- NULL
B <- 100
N <- 500
cutoff <- 1
models <- expand.grid(sigma = c(0.1, 0.5, 1),
                      x_dist = c("gauss", "exp", "unif", "gaussmix"),
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

    methods <- c("tsgauss", "emgg", "emgl")
    cat(": ")
    for (k in seq_along(methods))
    {
      cat(methods[k], ".. ")
      if (methods[k] == "tsgauss") {
        o <- tsgauss(dat$d, dat$w, cutoff)
      } else if (methods[k] == "emgg") {
        o <- emparam(dat$d, dat$w, cutoff,
                     x_dist = "gauss", u_dist = "gauss", verbose = FALSE)
      } else if (methods[k] == "emgl") {
        o <- emparam(dat$d, dat$w, cutoff,
                     x_dist = "gauss", u_dist = "lap", verbose = FALSE)
      }
      tmp <- data.frame(
        data_id = i, sigma = sigma, x_dist = x_dist,
        method = methods[k],
        param = names(o$estimate), estimate = o$estimate, stderr = o$stderr,
        convergence = o$convergence,
        stringsAsFactors = FALSE)
      results <- rbind(results, tmp)
    }
  }
  cat("\n")
}


# ggplot(results, aes(factor(data_id), estimate, color = method)) +
#   theme_bw() + facet_wrap(~param) +
#   geom_boxplot()

ggplot(subset(results, param == "sigma" & data_id != 3 & data_id != 11),
       aes(factor(data_id), estimate, color = method)) +
  theme_bw() +
  geom_boxplot()
ggsave("examples/sim-box.pdf", width = 9, height = 3)

o <- group_by(results, param, data_id, method) %>%
  summarize(sigma = mean(sigma),
            mean_est = mean(estimate),
            sd_est = sd(estimate),
            mean_se = mean(stderr)) %>%
  as.data.frame()
print(o)
ggplot(subset(o, param == "sigma"), aes(sd_est, mean_se, color = method)) +
  facet_wrap(~ sigma) +
  geom_point()

