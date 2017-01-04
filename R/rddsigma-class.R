#' @export
summary.rddsigma <- function(object, ...)
{
  if (object$model == "tsgauss") {
    model <- "Two Step Gaussian"
    x_dist <- "Gaussian"
    u_dist <- "Gaussian"
  } else if (object$model == "emparam") {
    model <- "EM Parametric"
    if (object$x_dist == "gauss") {
      x_dist <- "Gaussian"
    } else {
      x_dist <- object$x_dist
      warning("unsupported x_dist. can be fake 'rddsigma'")
    }
    if (object$u_dist == "gauss") {
      u_dist <- "Gaussian"
    } else if (object$u_dist == "lap") {
      u_dist <- "Laplace"
    } else {
      u_dist <- object$u_dist
      warning("unsupported u_dist. can be fake 'rddsigma'")
    }
  }

  est <- object$estimate
  se <- object$stderr
  n <- object$nobs
  zval <- est / se
  p <- pnorm(abs(zval), lower.tail = FALSE)
  coef_table <- cbind(est, se, zval, p)
  colnames(coef_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|t|)")
  rownames(coef_table) <- names(est)
  cat("* RDD sigma Estimate *\n\n")
  printCoefmat(coef_table)
  cat("\n")
  cat(" n obs  : ", n, "\n")
  cat(" Method : ", model, "\n")
  cat(" x dist : ", x_dist, "\n")
  cat(" u dist : ", u_dist, "\n")
  cat("\n")
}

#' @export
print.rddsigma <- function(x, ...) { summary(x) }

