library(testthat)
library(rddsigma)

context("Valid input")

test_that("valid input", {
  dat <- gen_data(500, 0.2, 0)

  expect_error(estimate_sigma(dat$d, dat$w, 0, method="foo"))
  expect_error(estimate_sigma(dat$d, dat$w, 0, method=pi))

})