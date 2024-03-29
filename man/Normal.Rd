% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normal.R
\name{Normal}
\alias{Normal}
\title{Normal Estimate}
\usage{
Normal(data, cL)
}
\arguments{
\item{data}{Data set/sample}

\item{cL}{confidence interval}
}
\value{
estimated mean, estimated standard deviation & confidence interval
}
\description{
Estimation of Mean, Standard Deviation and confidence interval
of Normal Distribution
}
\examples{
Normal(rnorm(2000,25,4),.95)
}
