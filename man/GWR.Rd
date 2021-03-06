% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GWR.R
\name{GWR}
\alias{GWR}
\title{GWR
to be documented}
\usage{
GWR(Y,XV,XX,YY,S,H, kernels, type = "GD", minv = 1, maxknn = 500,
 NmaxDist = 6000,SE=FALSE, isgcv, TIME, decay,interceptv=TRUE,W=NULL,
 betacor=FALSE,remove_local_outlier=FALSE,outv=0,doMC=FALSE,ncore=1,
 Wh=NULL,xratiomin=10e-10)
}
\arguments{
\item{Y}{A vector of response}

\item{XV}{A matrix with covariates with stationnary parameters}

\item{XX}{A matrix with all covariates (XC,XV)}

\item{YY}{A vector of first stage Y for some models (see MGWRSAR code)}

\item{S}{A matrix with variables used in kernel}

\item{H}{A vector of bandwidths}

\item{kernels}{A vector of kernel types}

\item{type}{Type of Genelarized kernel product ('GD' only spatial,'GDC' spatial + a categorical variable,'GDX' spatial + a continuous variable,'GDT' spatial + a time index, and other combinations 'GDXXC','GDTX',...)}

\item{minv}{Minimum number of non null weight}

\item{maxknn}{If n >NmaxDist how many column with dense weight matrix (max number of neighbours)}

\item{NmaxDist}{Maximum number of observation for computing dense weight matrix}

\item{SE}{If standard error are computed, default FALSE}

\item{isgcv}{leave one out cross validation, default FALSE.}

\item{TIME}{Use rigth truncated kernel for time index kernel}

\item{decay}{time decay}

\item{interceptv}{Intercept spatially varying, default FALSE}

\item{W}{A weight matrix for spatial autocorrelation}

\item{betacor}{Do a tuncation of spatial autocorelation if absolute value larger than 1.}

\item{remove_local_outlier}{Remove local outlier}

\item{outv}{A treshold for removing local outlier}

\item{doMC}{doParallel parallelization}

\item{ncore}{Number of cores for parallelization}

\item{Wh}{A matrix of weights for local estimation}

\item{xratiomin}{A treshold parameters for removing obs with not enough positive weigths for local regression}
}
\value{
a list of object for MGWRSAR wrapper
}
\description{
GWR
to be documented
}
\keyword{internal}
