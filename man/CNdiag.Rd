% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elicit_functions.R
\name{CNdiag}
\alias{CNdiag}
\title{Function to check condition number diagnostic.}
\usage{
CNdiag(X)
}
\arguments{
\item{X}{Design matrix}
}
\value{
a scalar giving the condition number of the rescaled design matrix
}
\description{
This function calculates the condition number of  the rescaled \eqn{n x
p} design matrix \eqn{X} such that each column has unit length.
}
\examples{
X <- matrix(rnorm(16), nrow = 4)
CNdiag(X)
}
