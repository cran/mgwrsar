\name{mgwrsar-internal}
\alias{INST_C}
\alias{PhWY_C}
\alias{Proj_C}
\alias{QRcpp2_C}
\alias{Sl_C}
\alias{gwr_beta_pivotal_qrp_cpp}
\alias{gwr_beta_univar_cpp}
\alias{mgwr_beta_pivotal_qrp_mixed_cpp}
\alias{ApproxiW}
\alias{mod}
\alias{get_index_mahalanobis_dual_rcpp}

\title{Internal Functions for mgwrsar Package}
\description{
  Internal functions used by the mgwrsar package. These functions are exported for
  technical reasons but are not intended to be called directly by the typical user.
}
\usage{
INST_C(A, B, C, D, nthreads = 1L)
PhWY_C(A, B, C, D, nthreads = 1L)
Proj_C(A, B, nthreads = 1L)
QRcpp2_C(A, B, C, nthreads = 1L)
Sl_C(A, B, C, D, nthreads = 1L)
gwr_beta_pivotal_qrp_cpp(X, y, XV, indexG, Wd, TP, get_ts = FALSE,
                          get_s = FALSE, get_Rk = FALSE, get_se = FALSE,
                          nthreads = 1L)
gwr_beta_univar_cpp(y, x, XV, indexG, Wd, TP, get_ts = FALSE, get_s = FALSE,
                    get_se = FALSE, nthreads = 1L)
mgwr_beta_pivotal_qrp_mixed_cpp(XV, y, XC, indexG, Wd, TP, get_ts = FALSE,
                                get_s = FALSE, get_Rk = FALSE, get_se = FALSE,
                                nthreads = 1L)

ApproxiW(W, TP, n, nthreads = 1L)
mod(Y, X, W, X_s, Y_s, S, type, b2sls, get_ts, get_s, nthreads = 1L)
get_index_mahalanobis_dual_rcpp(DS, DTM, nthreads = 1L)
}
\arguments{
  \item{A, B, C, D, X, y, XV, XC, indexG, Wd, TP, get_ts, get_s, get_Rk, get_se, x, W, n, X_s, Y_s, S, type, b2sls, DS, DTM, nthreads}{Internal arguments.}
}
\keyword{internal}
