## post simulation study

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(reshape2)


filelist <- dir("examples/sim-files/", "^simres", full.names = TRUE)
x <- NULL
for (fn in filelist)
{
  cat(fn, "...\n")
  a <- read_csv(fn)
  a$stderr <- as.numeric(a$stderr)
  x <- rbind(x, a)
}

y <- filter(x, !grepl("decon", estimator))


x %>% filter(variable == "sigma") %>%
  group_by(estimator, setup_id, convergence) %>%
  summarize(n = length(estimate)) %>%
  spread(key = convergence, value = n, fill = 0) %>%
  as.data.frame()

with(x, unique(paste(setup_id, sigma, x_dist, u_dist)))


g1 <- ggplot(filter(y, variable == "sigma", convergence == 0, sigma == 0.2),
             aes(factor(setup_id), estimate, color = estimator)) +
  geom_boxplot()
g2 <- ggplot(filter(y, variable == "sigma", convergence == 0, sigma == 1.2),
             aes(factor(setup_id), estimate, color = estimator)) +
  geom_boxplot()
grid.arrange(g1, g2)




