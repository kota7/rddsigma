library(rddme)
library(ggplot2)
library(dplyr)


## todo: parallelize this?

set.seed(123)

results <- NULL
B <- 500
N <- 1000
cutoff <- 0
models <- expand.grid(sigma = c(0.3, 1),
                      x_dist = c("gauss", "gaussmix"),
                      stringsAsFactors = FALSE)

for (i in 1:nrow(models))
{
  sigma <- models$sigma[i]
  x_dist <- models$x_dist[i]
  cat(i, "/", nrow(models), ": sigma =", sigma, ", x_dist =", x_dist, "\n")
  for (b in 1:B)
  {
    cat(sprintf("\r  %4d/%4d", b, B))
    dat <- gen_data(N, sigma, cutoff)
    o1 <- tsgauss(dat$d, dat$w, cutoff)

    tmp <- data.frame(
      data_id = i, sigma = sigma, x_dist = x_dist, method = "2 step gauss",
      param = names(o1$estimate), estimate = o1$estimate, stderr = o1$stderr,
      stringsAsFactors = FALSE)
    results <- rbind(results, tmp)
  }
  cat("\n")
}
ggplot(subset(results, param == "sigma"),
       aes(factor(data_id), estimate, color = method)) +
  theme_bw() +
  geom_boxplot()

ggplot(results, aes(factor(data_id), estimate, color = method)) +
  theme_bw() + facet_wrap(~param) +
  geom_boxplot()

group_by(results, param, data_id, method) %>%
  summarize(sigma = mean(sigma), mean(estimate), sd(estimate), mean(stderr))

