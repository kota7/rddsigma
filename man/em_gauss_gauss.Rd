% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_gauss_gauss.R
\name{em_gauss_gauss}
\alias{em_gauss_gauss}
\title{EM Algorithm Estimator for Gaussian-Gaussian Model}
\usage{
em_gauss_gauss(d_vec, w_vec, cutoff, reltol = 1e-06, maxit = 200L,
  integrate_options = list(), quiet = FALSE, ...)
}
\arguments{
\item{d_vec}{binary integer vector of assignment}

\item{w_vec}{numeric vector of observed running variable}

\item{cutoff}{threshold value for assignment}

\item{reltol}{relative tolerance requied}

\item{maxit}{maximum number of iteration}

\item{integrate_options}{controls to pass to \code{\link{integrate}}}

\item{quiet}{if false, progress is reported}

\item{...}{currently not used}
}
\value{
List
}
\description{
old version. kept for debugging.
}
\examples{
\dontrun{
dat <- gen_data(500, 0.2, 0)
em_gauss_gauss(dat$d, dat$w, 0)
}
}
\keyword{internal}
