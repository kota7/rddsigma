library(rddme)
library(microbenchmark)
library(ggplot2)

dat <- gen_data(1e+3, 0.2, 0)

m <- microbenchmark(
  tsgauss_r = rddme:::tsgauss_r(dat$d, dat$w, 0),
  tsgauss_cpp = tsgauss(dat$d, dat$w, 0)
)
print(m)
ggplot(as.data.frame(m), aes(expr, time*1e-6)) +
  theme_bw() + scale_y_log10() + ylab("time in millseconds") +
  geom_boxplot()
