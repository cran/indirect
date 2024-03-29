% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elicit_functions.R
\name{elicitPt}
\alias{elicitPt}
\title{Function to create or update elicitation at a given design point.}
\usage{
elicitPt(
  Z,
  design.pt = NULL,
  lower.CI.bound = NA,
  median = NA,
  upper.CI.bound = NA,
  CI.prob = NULL,
  comment = " "
)
}
\arguments{
\item{Z}{list of \code{design} with entries: \code{theta}, a \eqn{n x 4} 
matrix with columns that give lower, median and upper quantiles of the 
central credible interval followed by the probability \code{CI.prob} 
allocated to the interval; \code{link}, the link function used; and 
\code{target}. This list object is created by \code{\link{designLink}}}

\item{design.pt}{single integer that denotes design point of interest}

\item{lower.CI.bound}{scalar that gives the lower bound of the central 
credible interval, default \code{NA}.}

\item{median}{scalar value, default \code{NA}}

\item{upper.CI.bound}{scalar that gives the upper bound of the central 
credible interval, default \code{NA}.}

\item{CI.prob}{numeric, a fraction between 0 and 1 that defines probability 
attributed to central credible interval. For example, 1/2 for quartiles or
1/3 for tertiles. Default \code{NULL} uses the initial \code{CI.prob} as
defined by \code{\link{designLink}}.}

\item{comment}{character, ASCII text providing contributed commentary associated 
with elicitation design point. It is recommended to avoid special characters
such as quotation marks etc.}
}
\value{
\code{Z}, a list of \code{design} with entries: \code{theta}, a 
  \eqn{n x 4} matrix with columns that give lower, median and upper quantiles
  followed by \code{CI.prob}  with updated entries for row specified by
  argument \code{design.pt}; \code{link}, the link function used; and
  \code{target}.
}
\description{
Function to create or update elicitation at a given design point.
}
\examples{
X <- matrix(c(1, 1, 0, 1), nrow = 2) # design
Z <- designLink(design = X)
Z <- elicitPt(Z, design.pt = 1,
  lower.CI.bound = -1,
  median = 0,
  upper.CI.bound = 1,
  comment = "A completed elicitation scenario.")
}
