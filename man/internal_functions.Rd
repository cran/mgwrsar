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
INST_C(A, B, C, D)
PhWY_C(A, B, C, D)
Proj_C(A, B)
QRcpp2_C(A, B, C)
Sl_C(A, B, C, D)
gwr_beta_pivotal_qrp_cpp(X, y, XV, indexG, Wd, TP, get_ts = FALSE, get_s = FALSE,
                         get_Rk = FALSE, get_se = FALSE)
gwr_beta_univar_cpp(y, x, XV, indexG, Wd, TP, get_ts = FALSE, get_s = FALSE)
mgwr_beta_pivotal_qrp_mixed_cpp(XV, y, XC, indexG, Wd, TP, get_ts = FALSE,
                                get_s = FALSE, get_Rk = FALSE, get_se = FALSE)

ApproxiW(W, TP, n)
mod(Y, X, W, X_s, Y_s, S, type, b2sls, get_ts, get_s)
get_index_mahalanobis_dual_rcpp(DS, DTM)
}
\arguments{
  \item{A, B, C, D, X, y, XV, XC, indexG, Wd, TP, get_ts, get_s, get_Rk, get_se, x, W, n, X_s, Y_s, S, type, b2sls, DS, DTM}{Internal arguments.}
}
\keyword{internal}
