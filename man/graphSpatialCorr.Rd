% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.r
\name{graphSpatialCorr}
\alias{graphSpatialCorr}
\title{Convert a factor to numeric}
\usage{
graphSpatialCorr(Y, w, time, RHOT = NA)
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
