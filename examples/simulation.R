library(rddsigma)
library(ggplot2)
library(dplyr)



## generate and save simulation data
set.seed(123)
setup <- expand.grid(sigma = c(0.2, 1.2),
                     x_dist = c("gauss", "exp"),
                     u_dist = c("gauss", "lap"),
                     B = 100L, N = 500L, cutoff = 1,
                     stringsAsFactors = FALSE)

out <- NULL
cat("Generating data...\n")
for (i in 1:nrow(setup))
{
  cat(sprintf("\r%4d/%4d...", i, nrow(setup)))
  x <- gen_data(n = setup$N[i] * setup$B[i],
                sigma = setup$sigma[i], cutoff = setup$cutoff[i],
                u_dist = setup$u_dist[i], x_dist = setup$x_dist[i]) %>%
    as.data.frame()
  sample_id <- lapply(1:setup$B[i], rep, setup$N[i]) %>% unlist()
  obs_id <- rep(1:setup$N[i], setup$B[i])
  x <- cbind(data.frame(setup_id = i, sample_id = sample_id, obs_id = obs_id),
             x)
  x <- mutate(x, sigma = setup$sigma[i], cutoff = setup$cutoff[i],
              u_dist = setup$u_dist[i], x_dist = setup$x_dist[i])
  out <- rbind(out, x)
}
cat(" done!\n")
path <- "examples/simulation-data.csv"
write.csv(out, path, row.names = FALSE)
cat("simulation data have been saved as:", path, "\n")




#
# for (i in 1:nrow(models))
# {
#   sigma <- models$sigma[i]
#   x_dist <- models$x_dist[i]
#   u_dist <- models$u_dist[i]
#   cat(i, "/", nrow(models), ": sigma =", sigma,
#       ", x_dist =", x_dist, ", u_dist =", u_dist, "\n")
#   for (b in 1:B)
#   {
#     cat(sprintf("\r  %4d/%4d", b, B))
#     dat <- gen_data(N, sigma, cutoff, u_dist = u_dist, x_dist = x_dist)
#
#     methods <- c("tsgauss", "emgg", "emgl", "emdecong", "emdeconl")
#     cat(": ")
#     for (k in seq_along(methods))
#     {
#       cat(methods[k], ".. ")
#       if (methods[k] == "tsgauss") {
#         o <- tsgauss(dat$d, dat$w, cutoff)
#       } else if (methods[k] == "emgg") {
#         o <- emparam(dat$d, dat$w, cutoff,
#                      x_dist = "gauss", u_dist = "gauss", verbose = FALSE)
#       } else if (methods[k] == "emgl") {
#         o <- emparam(dat$d, dat$w, cutoff,
#                      x_dist = "gauss", u_dist = "lap", verbose = FALSE)
#       } else if (methods[k] == "emdecong") {
#         o <- emdecon(dat$d, dat$w, cutoff, u_dist = "gauss", verbose = FALSE)
#       } else if (methods[k] == "emdeconl") {
#         o <- emdecon(dat$d, dat$w, cutoff, u_dist = "lap", verbose = FALSE)
#       }
#       tmp <- data.frame(
#         data_id = i, sigma = sigma, x_dist = x_dist,
#         method = methods[k],
#         param = names(o$estimate), estimate = o$estimate, stderr = o$stderr,
#         convergence = o$convergence,
#         stringsAsFactors = FALSE)
#       results <- rbind(results, tmp)
#     }
#   }
#   cat("\n")
# }
#
#
# # ggplot(results, aes(factor(data_id), estimate, color = method)) +
# #   theme_bw() + facet_wrap(~param) +
# #   geom_boxplot()
#
# ggplot(subset(results, param == "sigma" & data_id != 3 & data_id != 11),
#        aes(factor(data_id), estimate, color = method)) +
#   theme_bw() +
#   geom_boxplot()
# ggsave("examples/sim-box.pdf", width = 9, height = 3)
#
# o <- group_by(results, param, data_id, method) %>%
#   summarize(sigma = mean(sigma),
#             mean_est = mean(estimate),
#             sd_est = sd(estimate),
#             mean_se = mean(stderr)) %>%
#   as.data.frame()
# print(o)
# ggplot(subset(o, param == "sigma"), aes(sd_est, mean_se, color = method)) +
#   facet_wrap(~ sigma) +
#   geom_point()

