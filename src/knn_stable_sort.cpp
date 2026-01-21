#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>
#include <cmath>

Rcpp::List knn_stable_sort(const Rcpp::NumericMatrix& dist,
                           const Rcpp::IntegerMatrix& idx) {

  int nrow = dist.nrow();
  int ncol = dist.ncol();
  double factor = std::pow(10.0, 14);

  NumericMatrix dist_out(nrow, ncol);
  IntegerMatrix idx_out(nrow, ncol);

  for (int i = 0; i < nrow; i++) {

    // Create vector of pairs (rounded distance, index)
    std::vector<std::pair<double, int>> v(ncol);

    for (int j = 0; j < ncol; j++) {
      double d = dist(i, j);

      // Round to 14 digits to stabilize across platforms
      d = std::round(d * factor) / factor;

      v[j] = std::make_pair(d, idx(i, j));
    }

    // Stable sort: first by distance, then by index (tie-breaker)
    std::stable_sort(
      v.begin(), v.end(),
      [](const std::pair<double,int> &a,
         const std::pair<double,int> &b) {
        if (a.first < b.first) return true;
        if (a.first > b.first) return false;
        return a.second < b.second;  // tie-break by index
      }
    );

    // Write sorted results back
    for (int j = 0; j < ncol; j++) {
      dist_out(i, j) = v[j].first;
      idx_out(i, j)  = v[j].second;
    }
  }
  return List::create(
    _["dist"] = dist_out,
    _["idx"]  = idx_out
  );
}


Rcpp::List compute_DS_DT_cpp(
    const NumericMatrix& coords,
    const NumericVector& Time,
    const IntegerMatrix& indexG,
    bool cyclic,
    double cycling
) {
  int nQ = indexG.nrow();
  int k  = indexG.ncol();

  NumericMatrix DS_out(nQ, k);
  NumericMatrix DT_out(nQ, k);

  for (int ii = 0; ii < nQ; ii++) {

    int i = indexG(ii, 0) - 1;  // reference index (1-based → 0-based)

    double xi = coords(i, 0);
    double yi = coords(i, 1);
    double ti = Time[i];

    for (int j = 0; j < k; j++) {

      int jj = indexG(ii, j) - 1;

      // ---- spatial distance ----
      double dx = coords(jj, 0) - xi;
      double dy = coords(jj, 1) - yi;
      DS_out(ii, j) = std::sqrt(dx*dx + dy*dy);

      // ---- temporal distance ----
      double dt = std::fabs(Time[jj] - ti);

      if (cyclic) {
        dt = std::fmod(dt, cycling);
        dt = std::min(dt, cycling - dt);
      }

      DT_out(ii, j) = dt;
    }
  }

  return List::create(
    _["DS"] = DS_out,
    _["DT"] = DT_out
  );
}
