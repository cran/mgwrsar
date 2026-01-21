// ====================================================================
// Unified header for Rcpp / RcppArmadillo / RcppEigen coexistence
// ====================================================================
#pragma once

// Define a protection macro to avoid RcppArmadillo conflicts
#ifndef RCPP_ARMADILLO_H
#  include <RcppArmadillo.h>
#endif

#ifndef RCPP_EIGEN_H
#  include <RcppEigen.h>
#endif

#include <Rcpp.h>

using namespace Rcpp;
