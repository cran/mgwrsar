// ============================================================================
// local_var.cpp — Local weighted variances through the neighbour index
// ============================================================================
// Used by get_HKmin() in R/TDS_MGWR.R (control_tds$check_pairs).
//
// prep_w() returns its row-normalised weights as an n x NN matrix in neighbour
// order: Wd(i, j) is the weight of site indexG(i, j) (1-based) for site i. For
// each site i and each column k of X, this computes in one pass over the
// neighbours
//     s1_k = sum_j Wd(i, j) * X(indexG(i, j), k)
//     s2_k = sum_j Wd(i, j) * X(indexG(i, j), k)^2
// and returns the n x p matrix of local variances s2 - s1^2. No gather matrix
// is materialised. Only the first ncol(Wd) columns of indexG are read, so Wd
// may have fewer columns than indexG (truncated kernels).
//
// OpenMP is optional (MGWRSAR_OPENMP=1 at build time); without it the loop is
// sequential and `nthreads` is ignored.

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

Rcpp::NumericMatrix local_var_cpp(const Rcpp::NumericMatrix& Wd,
                                  const Rcpp::IntegerMatrix& indexG,
                                  const Rcpp::NumericMatrix& X,
                                  int nthreads) {
  const int n  = Wd.nrow();
  const int m  = Wd.ncol();
  const int p  = X.ncol();
  const int nx = X.nrow();
  const int gl = indexG.nrow();

  if (gl != n)
    Rcpp::stop("local_var_cpp: nrow(indexG) must equal nrow(Wd)");
  if (indexG.ncol() < m)
    Rcpp::stop("local_var_cpp: indexG must have at least ncol(Wd) columns");

  const double* w = Wd.begin();
  const int*    g = indexG.begin();
  const double* x = X.begin();

  // Index range checked up front: Rcpp::stop() must not be called from inside
  // a parallel region.
  for (int j = 0; j < m; ++j) {
    const int* gj = g + static_cast<std::size_t>(j) * gl;
    for (int i = 0; i < n; ++i) {
      const int idx = gj[i];
      if (idx < 1 || idx > nx)
        Rcpp::stop("local_var_cpp: indexG entries must lie in 1..nrow(X)");
    }
  }

  Rcpp::NumericMatrix out(n, p);
  double* o = out.begin();

  // X transposed (covariate index fastest): the p values of a gathered site
  // are then contiguous, instead of p reads one column apart.
  std::vector<double> xt(static_cast<std::size_t>(nx) * p);
  for (int k = 0; k < p; ++k)
    for (int i = 0; i < nx; ++i)
      xt[static_cast<std::size_t>(i) * p + k] = x[i + static_cast<std::size_t>(k) * nx];

  // Blocks of consecutive sites: within a block the neighbour columns of Wd
  // and indexG are read contiguously, and the block's accumulators (also
  // covariate-fastest) stay in cache. Blocks are independent, hence the
  // (optional) parallel loop.
  const int B = 256;
  const int nblocks = (n + B - 1) / B;
  if (nthreads < 1) nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for (int b = 0; b < nblocks; ++b) {
    const int i0 = b * B;
    const int i1 = std::min(n, i0 + B);
    const int nb = i1 - i0;
    std::vector<double> s1(static_cast<std::size_t>(nb) * p, 0.0);
    std::vector<double> s2(static_cast<std::size_t>(nb) * p, 0.0);
    for (int j = 0; j < m; ++j) {
      const double* wj = w + static_cast<std::size_t>(j) * n + i0;
      const int*    gj = g + static_cast<std::size_t>(j) * gl + i0;
      for (int i = 0; i < nb; ++i) {
        const double wij = wj[i];
        if (wij == 0.0) continue;
        const double* xi = xt.data() + static_cast<std::size_t>(gj[i] - 1) * p;
        double* a1 = s1.data() + static_cast<std::size_t>(i) * p;
        double* a2 = s2.data() + static_cast<std::size_t>(i) * p;
        for (int k = 0; k < p; ++k) {
          const double v = xi[k];
          a1[k] += wij * v;
          a2[k] += wij * v * v;
        }
      }
    }
    for (int i = 0; i < nb; ++i)
      for (int k = 0; k < p; ++k) {
        const double a = s1[static_cast<std::size_t>(i) * p + k];
        o[i0 + i + static_cast<std::size_t>(k) * n] =
          s2[static_cast<std::size_t>(i) * p + k] - a * a;
      }
  }
  return out;
}
