% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atds_gwr.R
\name{atds_gwr}
\alias{atds_gwr}
\title{atds_gwr Top-Down Scaling approach of GWR}
\usage{
atds_gwr(formula,data,coords,kernels='triangle',fixed_vars=NULL,
control_tds=list(nns=30),control=list(adaptive=TRUE,verbose=FALSE))
}
\arguments{
\item{formula}{a formula.}

\item{data}{a dataframe.}

\item{coords}{default NULL, a dataframe or a matrix with coordinates}

\item{kernels}{A vector containing the kernel types. Possible types:
triangle ("triangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane").}

\item{fixed_vars}{a vector with the names of spatiallay constant
coefficient for mixed model. All other variables present in formula
are supposed to be spatially varying. If empty or NULL (default),
all variables in formula are supposed to be spatially varying.}

\item{control_tds}{list of extra control arguments for tds_mgwr model - see tds_gwr Help}

\item{control}{list of extra control arguments for MGWRSAR wrapper - see MGWRSAR Help}
}
\description{
This function performs a Geographically Weighted Regression (GWR) using
a top-down scaling approach, adjusting GWR coefficients with a progressively
decreasing bandwidth as long as the AICc criterion improves.
}
\seealso{
tds_mgwr, gwr_multiscale, MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar.
}
