% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mash-evaluation.R
\name{scale_cov}
\alias{scale_cov}
\title{Scale each covariance matrix in list Ulist by a scalar in
vector grid}
\usage{
scale_cov(Ulist, grid)
}
\arguments{
\item{Ulist}{a list of matrices}

\item{grid}{a vector of scaling factors (standard deviaions)}
}
\value{
a list with length length(Ulist)*length(grid)
}
\description{
This is an internal (non-exported) function. This help
page provides additional documentation mainly intended for
developers and expert users.
}
\keyword{internal}
