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
path <- "examples/sim-files/simulation-data.csv"
write.csv(out, path, row.names = FALSE)
cat("simulation data have been saved as:", path, "\n")


