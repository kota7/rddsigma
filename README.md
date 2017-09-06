
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/kota7/rddsigma.svg?branch=master)](https://travis-ci.org/kota7/rddsigma)

rddsigma
========

Sharp Regression Discontinuity Design with Running Variable Measured with Error

This package provides estimators for the standard deviation of measurement errors in running variables for sharp regression discontinuity designs. The estimation requires the data for observed (mismeasured) running variables and assigment variables.

Installation
------------

Install from GitHub using `devtools::install_github` function.

``` r
devtools::install_github("kota7/rddsigma")
```

Quick Start
-----------

Use `estimate_sigma` function.

``` r
library(rddsigma)
set.seed(875)

# generate simulation data
dat <- gen_data(500, 0.3, 0)

estimate_sigma(dat$d, dat$w, cutoff=0)
#> * RDD sigma Estimate *
#> 
#>        Estimate Std. Error z value Pr(>|t|)    
#> sigma  0.286483   0.025592  11.194   <2e-16 ***
#> mu_x  -0.010614   0.046339  -0.229   0.4094    
#> sd_x   0.995779   0.034677  28.716   <2e-16 ***
#> sd_w   1.036170   0.032767  31.623   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  n obs  :  500 
#>  Method :  Two Step Gaussian 
#>  x dist :  Gaussian 
#>  u dist :  Gaussian 
#>  value  :  -0.195047 
#>  convergence: yes
```

Estimation Methods
------------------

The package support two estimators for estimating `sigma`. The first estimator assumes that both the true runnning variable and the measurement error follow the Gaussian distribution. This is the default estimator of `estimate_sigma` function, or specified explicitly by setting `method="tsgauss"`

``` r
estimate_sigma(dat$d, dat$w, cutoff=0, method="tsgauss")
#> * RDD sigma Estimate *
#> 
#>        Estimate Std. Error z value Pr(>|t|)    
#> sigma  0.286483   0.025592  11.194   <2e-16 ***
#> mu_x  -0.010614   0.046339  -0.229   0.4094    
#> sd_x   0.995779   0.034677  28.716   <2e-16 ***
#> sd_w   1.036170   0.032767  31.623   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  n obs  :  500 
#>  Method :  Two Step Gaussian 
#>  x dist :  Gaussian 
#>  u dist :  Gaussian 
#>  value  :  -0.195047 
#>  convergence: yes
```

The second estimator relaxes the Gaussian assumption and allows the true runnning variable and the measurement error to follow some parametric distribution. Currently, the package supports the Gaussian distribution for the running variable and the Gaussian and Laplace distributions for the measurement error. This estimator is called with `method="emparam"`

``` r
estimate_sigma(dat$d, dat$w, cutoff=0, method="emparam", x_dist="gauss", u_dist="gauss")
#> * RDD sigma Estimate *
#> 
#>        Estimate Std. Error z value Pr(>|t|)    
#> sigma  0.281109   0.024359  11.540   <2e-16 ***
#> mu_x  -0.010614   0.045546  -0.233   0.4079    
#> sd_x   0.998358   0.031221  31.978   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  n obs  :  500 
#>  Method :  EM Parametric 
#>  x dist :  Gaussian 
#>  u dist :  Gaussian 
#>  value  :  -1.650564 
#>  convergence: yes
estimate_sigma(dat$d, dat$w, cutoff=0, method="emparam", x_dist="gauss", u_dist="lap")
#> * RDD sigma Estimate *
#> 
#>        Estimate Std. Error z value  Pr(>|t|)    
#> sigma  0.299957   0.038387  7.8141 2.769e-15 ***
#> mu_x  -0.010614   0.045546 -0.2330    0.4079    
#> sd_x   0.991398   0.034457 28.7724 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  n obs  :  500 
#>  Method :  EM Parametric 
#>  x dist :  Gaussian 
#>  u dist :  Laplace 
#>  value  :  -1.653077 
#>  convergence: yes
```
