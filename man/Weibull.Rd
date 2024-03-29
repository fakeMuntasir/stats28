% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Weibull.R
\name{Weibull}
\alias{Weibull}
\title{Title Weibull distibution parameters estimate.}
\usage{
Weibull(data, cL)
}
\arguments{
\item{data}{data set / sample ( as a list ) to estimate the parameters and confidence interval.}

\item{cL}{confidence level.}
}
\value{
Point estimates the shape and scale parameter, Confidence interval for the shape and scale parameter.
}
\description{
Title Weibull distibution parameters estimate.
}
\examples{
(rweibull(0.5,.5),.95)
(rweibull(10,2),.90)
}
