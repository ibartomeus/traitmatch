% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/useful_functions.R
\name{Useful_functions}
\alias{Useful_functions}
\alias{weighted_mean}
\title{weighted mean and sd}
\usage{
weighted_mean(x, w)
}
\arguments{
\item{x}{vector to compute the mean over.}

\item{w}{a vector of same length with the weigths of each element.}
}
\value{
a value with the weighted mean or sd.
}
\description{
just weighted mean and sd.
}
\examples{
weighted_mean(c(2,3,4), c(1,1,2))
weighted_sd(c(2,3,4), c(2,1,2))

}

