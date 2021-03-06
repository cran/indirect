% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting_functions.R
\name{checkX}
\alias{checkX}
\title{Helper function that checks for sensible covariate matrix.}
\usage{
checkX(X)
}
\arguments{
\item{X}{numeric matrix of covariates, \eqn{n} design points by \eqn{p}
covariates, for a given model and design points.}
}
\value{
throws an error if not full rank.
}
\description{
Helper function that checks for sensible covariate matrix.
}
