% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdsr.r
\name{loglikTVRhoVarTrd}
\alias{loglikTVRhoVarTrd}
\title{Convert a factor to numeric}
\usage{
loglikTVRhoVarTrd(
  Y,
  w,
  omegaRho = 0,
  aRho = 0.01,
  bRho = 0.8,
  f1Rho = atanh(0.4),
  omegaVar = 0,
  aVar = 0.01,
  bVar = 0.8,
  f1Var = log(1),
  omegaTrd = 0,
  aTrd = 0.01,
  bTrd = 0.8,
  f1Trd = 0,
  density = "normal",
  result = "loglik"
)
}
\arguments{
\item{x}{A vector containing a factor with numeric levels}
}
\value{
The input factor made a numeric vector
}
\description{
Convert a factor with numeric levels to a non-factor
}
\examples{
x <- factor(c(3, 4, 9, 4, 9), levels=c(3,4,9))
fac2num(x)

}
