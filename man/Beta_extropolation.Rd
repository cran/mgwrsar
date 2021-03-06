% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shepard_smooth_16tg.R
\name{Beta_extropolation}
\alias{Beta_extropolation}
\title{An extrapolation method using shepard smooth with a given number of neighbors.}
\usage{
Beta_extropolation(coord=coord,k=16,O=NULL,type='SxS')
}
\arguments{
\item{coord}{a matrix of coordinates}

\item{k}{The number of neighbors for extrapolation}

\item{O}{An index vector for outsample}

\item{type}{Insample/outsample size configuration : 'S+OxS+O' or 'SxS'(default)}
}
\value{
A weight dgCMatrix of size $type$
}
\description{
An extrapolation method using shepard smooth with a given number of neighbors.
}
\keyword{internal}
