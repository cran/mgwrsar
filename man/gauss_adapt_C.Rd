% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gauss_adapt_C.R
\name{gauss_adapt_C}
\alias{gauss_adapt_C}
\title{Adaptive gaussian kernel, RcppEigen version}
\usage{
gauss_adapt_C(d, h)
}
\arguments{
\item{d}{a vector of distance}

\item{h}{bandwidth size expressed in number of neighbors}
}
\value{
a vector of weight
}
\description{
Adaptive gaussian kernel, RcppEigen version
}
\examples{
\donttest{
w=gauss_adapt_C(-30:30,20)
plot(-30:30,w,type='l')
abline(v=-10)
abline(v=10)
}
}
