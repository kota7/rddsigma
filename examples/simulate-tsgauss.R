library(rddsigma)
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)

## Simulation of Two step gauss estimator


## load data
x <- read_csv("examples/sim-files/simulation-data.csv")


## do it!
out <- NULL
cat("Simulating for Two-step Gauss Estimator...\n")
for (i in unique(x$setup_id))
{
  y <- filter(x, setup_id == i)
  for (b in unique(y$sample_id))
  {
    cat(sprintf("\rsetup = %2d/%2d, sample = %4d/%4d ...",
                i, length(unique(x$setup_id)),
                b, length(unique(y$sample_id))))
    z <- filter(y, sample_id == b)
    u <- select(z, setup_id, sample_id, sigma, cutoff, u_dist, x_dist) %>%
      distinct()
    ## data generation parameters must be unique within a sample
    stopifnot(nrow(u) == 1L)

    ## estimate
    fit <- tsgauss(z$d, z$w, z$cutoff[1])

    ## stack
    o <- data.frame(variable = names(fit$estimate),
                    estimate = fit$estimate,
                    stderr = fit$stderr,
                    convergence = fit$convergence)
    o <- cbind(data.frame(estimator = "tsgauss", stringsAsFactors = FALSE), o)
    o <- merge(o, u)

    out <- rbind(out, o)
  }
}

## save it!
write.csv(out, "examples/sim-files/simres-tsgauss.csv", row.names = FALSE)


## visualize it!
g1 <- ggplot(filter(out, sigma == 0.2, variable == "sigma"),
             aes(factor(setup_id), estimate)) + geom_boxplot()
g2 <- ggplot(filter(out, sigma == 1.2, variable == "sigma"),
             aes(factor(setup_id), estimate)) + geom_boxplot()
grid.arrange(g1, g2)
