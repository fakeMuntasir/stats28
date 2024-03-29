% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Poisson.R
\name{Poisson}
\alias{Poisson}
\title{Poisson Estimate}
\usage{
Poisson(data, cL)
}
\arguments{
\item{data}{Data set/sample}

\item{cL}{confidence interval}
}
\value{
estimated parameter lambda, confidence interval
}
\description{
Estimation of the parameter & confidence interval from the sample
of a Poisson Distribution
}
\examples{
poisson(rpois(100,.5),.95)
poisson(rpois(1000,.45),.95)
}
