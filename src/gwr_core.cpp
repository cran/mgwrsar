// ---------------------------------------------------------------------------
// gwr_core.cpp — GWR / MGWR Module for mgwrsar
// ---------------------------------------------------------------------------
// This file contains C++ (RcppArmadillo) functions for computing
// locally weighted coefficients:
//  - gwr_beta_univar_cpp: fast univariate version
//  - gwr_beta_pivotal_qrp_cpp: multivariate version using pivoted QR decomposition
//  - get_index_mahalanobis_dual_rcpp: combined space-time distances
// ---------------------------------------------------------------------------

#include "RcppIncludesArmadillo.h"
#include <fstream>
#include <iomanip>

// ---------------------------------------------------------------------------
// UTILITY FUNCTION — Effective rank of a matrix R (from QR decomposition)
// ---------------------------------------------------------------------------
static inline arma::uword eff_rank_from_R(const arma::mat& R) {
  if (R.n_rows == 0 || R.n_cols == 0) return 0;
  arma::vec d = arma::abs(R.diag());
  const double tol = std::max(R.n_rows, R.n_cols) * d.max() *
    std::numeric_limits<double>::epsilon();
  arma::uword r = 0;
  for (; r < d.n_elem; ++r)
    if (d[r] <= tol) break;
  return r;
}

// Apply Q' = H_r ... H_1 (the first r Householder reflectors of a LAPACK QR
// factorisation held in A and tau: v_k(k) = 1, v_k(i > k) = A(i, k)) to x.
// Used for Q'y and for the focal row of Q (Q' e_f) without forming Q.
static inline void apply_qt(const arma::mat& A, const arma::vec& tau,
                            const arma::uword r, arma::vec& x) {
  const arma::uword m = A.n_rows;
  for (arma::uword k = 0; k < r; ++k) {
    const double* vk = A.colptr(k);
    double s = x[k];
    for (arma::uword i = k + 1; i < m; ++i) s += vk[i] * x[i];
    s *= tau[k];
    x[k] -= s;
    for (arma::uword i = k + 1; i < m; ++i) x[i] -= s * vk[i];
  }
}

// ---------------------------------------------------------------------------
// Neighbour layout and stderr silencing shared by the GWR cores
// ---------------------------------------------------------------------------
// The wrappers pass indexG and Wd transposed (NN x nTP), so that the neighbour
// list of focal point z is the contiguous column z. Reading rows of the
// column-major n x NN matrices was cache-hostile and dominated the run time.

// Redirects Rcpp::Rcerr to a sink for the lifetime of the object (RAII), so
// the stream is restored even when the guarded solve throws.
struct RcerrSilencer {
  std::streambuf* old;
  explicit RcerrSilencer(std::streambuf* sink) : old(Rcpp::Rcerr.rdbuf(sink)) {}
  ~RcerrSilencer() { Rcpp::Rcerr.rdbuf(old); }
};

// ---------------------------------------------------------------------------
// Internal C++ function: gwr_beta_univar_core()
// ---------------------------------------------------------------------------
// Single fused pass over the neighbours of each focal point: no gather copies,
// no per-point heap allocation.
Rcpp::List gwr_beta_univar_core(const arma::vec& y,
                                const arma::vec& x,
                                const arma::mat& XV,
                                const arma::Mat<int>& idxT,
                                const arma::mat& WdT,
                                const arma::uvec& TP,
                                const bool get_ts,
                                const bool get_s,
                                const bool get_se) {

  const arma::uword nTP = TP.n_elem;
  const arma::uword n   = x.n_elem;
  const arma::uword NN  = idxT.n_rows;

  arma::vec Betav(nTP, arma::fill::zeros);
  arma::vec SEV(nTP, arma::fill::zeros);
  arma::vec TS(nTP, arma::fill::zeros);
  arma::mat Shat;
  if (get_s) Shat = arma::zeros(nTP, n);

  std::vector<arma::uword> rows(NN);
  std::vector<double> w_loc(NN);

  for (arma::uword z = 0; z < nTP; ++z) {
    const int*    id = idxT.colptr(z);
    const double* wz = WdT.colptr(z);

    arma::uword m = 0;
    double Sxx = 0.0, Sxy = 0.0;
    for (arma::uword j = 0; j < NN; ++j) {
      const double wj = wz[j];
      const int    k1 = id[j];               // 1-based; NA_INTEGER < 1
      if (!(wj > 1e-12) || k1 < 1 || arma::uword(k1) > n) continue;
      const arma::uword k = arma::uword(k1) - 1u;
      const double wx = wj * x[k];
      Sxx += wx * x[k];
      Sxy += wx * y[k];
      rows[m]  = k;
      w_loc[m] = wj;
      ++m;
    }

    if (m == 0) continue;

    const double beta = (Sxx == 0.0) ? 0.0 : Sxy / Sxx;
    Betav[z] = beta;

    if (get_se && Sxx > 0.0) {
      double rss_loc = 0.0;
      for (arma::uword j = 0; j < m; ++j) {
        const double res = y[rows[j]] - beta * x[rows[j]];
        rss_loc += (w_loc[j] * res) * res;
      }
      const double sigma2_loc = rss_loc / std::max(1.0, double(m) - 1.0);
      SEV[z] = std::sqrt(sigma2_loc / Sxx);
    }

    if ((get_ts || get_s) && Sxx != 0.0) {
      const double c = XV(TP[z] - 1u, 0) / Sxx;

      if (get_ts) {
        const arma::uword focal_id = TP[z] - 1u;
        for (arma::uword j = 0; j < m; ++j) {
          if (rows[j] == focal_id) {
            TS[z] = w_loc[j] * (x[rows[j]] * c);
            break;
          }
        }
      }

      if (get_s)
        for (arma::uword j = 0; j < m; ++j)
          Shat(z, rows[j]) = w_loc[j] * (x[rows[j]] * c);
    }
  }

  Rcpp::List out;
  out["Betav"] = Betav;
  if (get_se) out["SEV"] = SEV;
  if (get_ts) out["TS"] = TS;
  if (get_s)  out["Shat"] = Shat;
  return out;
}

// ---------------------------------------------------------------------------
// Updated Rcpp wrapper
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List gwr_beta_univar_cpp(const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector& x,
                               const Rcpp::NumericMatrix& XV,
                               const Rcpp::IntegerMatrix& indexG,
                               const Rcpp::NumericMatrix& Wd,
                               const Rcpp::IntegerVector& TP,
                               bool get_ts,
                               bool get_s,
                               bool get_se) { // Added argument

  if (Wd.ncol() < indexG.ncol() || Wd.nrow() < TP.size() || indexG.nrow() < TP.size())
    Rcpp::stop("indexG and Wd must have one row per target point and matching columns.");

  arma::vec ay(const_cast<double*>(y.begin()), y.size(), false);
  arma::vec ax(const_cast<double*>(x.begin()), x.size(), false);
  arma::mat aXV(const_cast<double*>(XV.begin()), XV.nrow(), XV.ncol(), false);
  arma::Mat<int> idxT = arma::Mat<int>(const_cast<int*>(indexG.begin()),
                                       indexG.nrow(), indexG.ncol(), false).t();
  arma::mat WdT = arma::mat(const_cast<double*>(Wd.begin()),
                            Wd.nrow(), Wd.ncol(), false).t();
  arma::uvec aTP = Rcpp::as<arma::uvec>(TP);

  return gwr_beta_univar_core(ay, ax, aXV, idxT, WdT, aTP, get_ts, get_s, get_se);
}
// ---------------------------------------------------------------------------
// Main function: gwr_beta_pivotal_qrp_core (multivariate pivoted QR)
// ---------------------------------------------------------------------------

// Xt is X transposed (p x n): the covariates of one observation are contiguous.
// Buffers are allocated once outside the loop. Each local system is solved by
// a column-pivoted Householder QR (LAPACK geqp3): the diagonal of R is then
// non-increasing, so the rank is read reliably and the non-estimable columns
// are the pivoted ones beyond the rank, whatever their position in the
// formula. Hat quantities use Q only: with Xw_r = Q_r Rr and x_f the focal
// row, S(f, j) = (sw_j / sw_f) q_f' q_j and s_ff = ||q_f||^2, which stays in
// [0, 1] even when the local system is nearly singular (R^{-1} is used only
// for the coefficients, their standard errors and the per-coefficient
// operator Rk). Without Shat/Rk, Q is never formed: Q'y and q_f come from the
// reflectors.
Rcpp::List gwr_beta_pivotal_qrp_core(
    const arma::mat& Xt,
    const arma::vec& y,
    const arma::mat& XV,
    const arma::Mat<int>& idxT,
    const arma::mat& WdT,
    const arma::uvec& TP,
    const bool get_ts,
    const bool get_s,
    const bool get_Rk,
    const bool get_se
) {
  using namespace arma;

  const uword nTP = TP.n_elem;
  const uword n   = Xt.n_cols;
  const uword p   = Xt.n_rows;
  const uword NN  = idxT.n_rows;

  mat Betav(nTP, p, fill::zeros);

  vec TS(nTP, fill::zeros);
  mat Shat; if (get_s)  Shat = mat(nTP, n, fill::zeros);
  cube Rk;  if (get_Rk) Rk   = cube(nTP, n, p, fill::zeros);

  mat SEV;  if (get_se) SEV  = mat(nTP, p, fill::zeros);

  std::vector<uword> rows(NN);
  std::vector<double> sw(NN);
  mat Xw, R, A;
  vec yw;
  std::stringstream sink;

  // Q_r is formed (orgqr) only when the hat rows (Shat) or Rk are requested
  const bool need_Q = get_s || get_Rk;
  vec tau(std::max<uword>(p, 1u));
  vec work(std::max<uword>(64u * p, 4u));   // geqp3 needs >= 3p + 1, orgqr >= r
  std::vector<blas_int> jpvt(std::max<uword>(p, 1u));
  std::vector<uword> perm(std::max<uword>(p, 1u));

  // ---- loop over focal points
  for (uword z = 0; z < nTP; ++z) {

    const int*    id = idxT.colptr(z);
    const double* wz = WdT.colptr(z);

    // keep only valid neighbors (and keep their weights aligned)
    uword n_loc = 0;
    for (uword j = 0; j < NN; ++j) {
      const int    k1 = id[j];               // 1-based; NA_INTEGER < 1
      const double w  = wz[j];
      if (k1 >= 1 && uword(k1) <= n && w > 0) {
        rows[n_loc] = uword(k1) - 1u;
        sw[n_loc]   = std::sqrt(w);
        ++n_loc;
      }
    }

    if (n_loc <= p) continue;

    Xw.set_size(n_loc, p);
    yw.set_size(n_loc);
    for (uword j = 0; j < n_loc; ++j) {
      const double* xk = Xt.colptr(rows[j]);
      for (uword c = 0; c < p; ++c) Xw.at(j, c) = xk[c] * sw[j];
      yw[j] = y[rows[j]] * sw[j];
    }

    // ---- column-pivoted QR of Xw: Xw P = Q R
    A = Xw;
    blas_int m_i = blas_int(n_loc), p_i = blas_int(p), lda = blas_int(n_loc);
    blas_int lwork_i = blas_int(work.n_elem), info = 0;
    for (uword c = 0; c < p; ++c) jpvt[c] = 0;
    lapack::geqp3(&m_i, &p_i, A.memptr(), &lda, jpvt.data(), tau.memptr(),
                  work.memptr(), &lwork_i, &info);
    if (info != 0) continue;
    R = trimatu(A.rows(0, p - 1));
    const uword r = eff_rank_from_R(R);
    if (r == 0u) continue;
    const mat Rr = R.submat(0, 0, r - 1, r - 1);
    for (uword c = 0; c < p; ++c) perm[c] = uword(jpvt[c]) - 1u;

    // ---- coefficients: Rr beta_r = Q_r' y, zeros on the non-estimable columns
    vec Qt_y = yw;
    apply_qt(A, tau, r, Qt_y);
    vec beta_r;
    try {
      RcerrSilencer quiet(sink.rdbuf());
      beta_r = solve(trimatu(Rr), Qt_y.head(r), solve_opts::fast);
    } catch (...) {
      continue;
    }
    vec beta(p, fill::zeros);
    for (uword k = 0; k < r; ++k) beta[perm[k]] = beta_r[k];
    Betav.row(z) = beta.t();

    // ---- local standard errors (only if requested)
    if (get_se) {
      vec sev_full(p, fill::zeros);
      try {
        vec yhat_w(n_loc, fill::zeros);
        for (uword k = 0; k < r; ++k) yhat_w += Xw.col(perm[k]) * beta_r[k];
        const vec resid_w = yw - yhat_w;
        const double sigma2 = dot(resid_w, resid_w) /
          std::max(1.0, double(n_loc) - double(r));
        const mat Rinverse = inv(trimatu(Rr));            // r x r
        const vec sev_loc = sqrt(sigma2 * sum(square(Rinverse), 1)); // diag(Rinv Rinv')
        for (uword k = 0; k < r; ++k) {
          const double v = sev_loc[k];
          sev_full[perm[k]] = (std::isfinite(v) ? v : 0.0);
        }
      } catch (...) {
        // keep zeros (legacy behavior)
      }
      SEV.row(z) = sev_full.t();
    }

    // ---- hat quantities at the focal point
    const uword focal_id = TP[z] - 1u;
    uword jf = n_loc;
    for (uword jj = 0; jj < n_loc; ++jj)
      if (rows[jj] == focal_id) { jf = jj; break; }
    if (jf == n_loc) continue;   // the focal point is not among its neighbours

    if (!need_Q) {
      if (get_ts) {
        // s_ff = ||q_f||^2 with q_f = (Q' e_f)(1:r), from the reflectors
        vec qf(n_loc, fill::zeros);
        qf[jf] = 1.0;
        apply_qt(A, tau, r, qf);
        TS[z] = dot(qf.head(r), qf.head(r));
      }
      continue;
    }

    // Q_r formed in place in the first r columns of A
    blas_int r_i = blas_int(r);
    lapack::orgqr(&m_i, &r_i, &r_i, A.memptr(), &lda, tau.memptr(),
                  work.memptr(), &lwork_i, &info);
    if (info != 0) continue;
    const mat Qr = A.cols(0, r - 1);
    const vec qf = Qr.row(jf).t();
    vec s_local = Qr * qf;                       // q_f' q_j for every neighbour
    for (uword j = 0; j < n_loc; ++j) s_local[j] *= sw[j] / sw[jf];

    if (get_ts) TS[z] = s_local[jf];
    if (get_s)
      for (uword j = 0; j < n_loc; ++j) Shat(z, rows[j]) = s_local[j];

    if (get_Rk) {
      // per-coefficient operator: B = P [Rr^{-1} Q_r'; 0], Rk(f, j, nx) = x_f[nx] B(nx, j) sw_j
      mat Br;
      try {
        RcerrSilencer quiet(sink.rdbuf());
        Br = solve(trimatu(Rr), Qr.t(), solve_opts::fast);   // r x n_loc
      } catch (...) {
        continue;
      }
      const rowvec x0 = XV.row(focal_id);
      const uword focal = focal_id;
      for (uword k = 0; k < r; ++k) {
        const uword nx = perm[k];
        const double x0nx = x0[nx];
        if (std::isnan(x0nx)) continue;
        for (uword j = 0; j < n_loc; ++j)
          Rk(focal, rows[j], nx) = x0nx * Br(k, j) * sw[j];
      }
    }
  }

  // ---- effective parameters
  double tS  = sum(TS);
  double edf = double(n) - tS;

  // ---- output list (legacy-compatible)
  Rcpp::List out;
  out["Betav"] = Betav;

  if (get_se) out["SEV"] = SEV; else out["SEV"] = R_NilValue;

  out["edf"] = edf;
  out["tS"]  = tS;

  if (get_ts) out["TS"] = TS;
  if (get_s)  out["Shat"] = Shat;
  if (get_Rk) out["Rk"] = Rk;

  return out;
}

// ============================================================================
// MGWR mixed core  (same API names as before)
// ============================================================================
// ============================================================================
// MGWR mixed core (variable part local, fixed part global on ZXc = (I - Sv)XC)
// Adds:
//  - returns Shat (variable-part only) when get_s
//  - returns TS (diag(Sv) at focal) when get_ts
//  - computes diag(Hfix) for fixed part on ZXc
//  - computes approx diag(S_total) = diag(Sv) + diag(Hfix) * (1 - diag(Sv))
//  - returns tS_var and tS_tot_approx (and keeps out["tS"] backward compatible)
//  - uses tS_tot_approx (recommended) in sigma2_diag denominator
// Layout as in gwr_beta_pivotal_qrp_core: XVt, XCt (covariates x n) and the
// transposed neighbour arrays idxT, WdT (NN x nTP); buffers allocated once.
// get_Rk only triggers the TS computation (no Rk is returned by this core).
// ============================================================================

Rcpp::List mgwr_beta_pivotal_qrp_mixed_core_new(
    const arma::mat& XV,
    const arma::mat& XVt,
    const arma::vec& y,
    const arma::mat& XC,
    const arma::mat& XCt,
    const arma::Mat<int>& idxT,
    const arma::mat& WdT,
    const arma::uvec& TP,
    const bool get_ts,
    const bool get_s,
    const bool get_Rk,
    const bool get_se
) {

  arma::uword nTP = TP.n_elem,
    n   = XV.n_rows,
    kv  = XV.n_cols,
    kc  = XC.n_cols,
    NN  = idxT.n_rows;

  arma::mat SY(nTP, kv, arma::fill::zeros);
  // local mapping of XC through the variable fit, one kv x kc slice per z
  arma::cube XCw_cube(kv, kc, nTP, arma::fill::zeros);

  // TS stores diag(Sv) at focal point (i=TP[z]-1) for each z
  arma::vec TS(nTP, arma::fill::zeros);

  arma::mat SEV(nTP, kv, arma::fill::zeros);

  arma::mat Shat;
  if (get_s) Shat = arma::zeros(nTP, n);

  std::vector<arma::uword> rows(NN);
  std::vector<double> sw(NN);
  arma::mat Xw, XCw, R, A;
  arma::vec yw;
  std::stringstream sink;
  arma::vec tau(std::max<arma::uword>(kv, 1u));
  arma::vec work(std::max<arma::uword>(64u * kv, 4u));
  std::vector<blas_int> jpvt(std::max<arma::uword>(kv, 1u));
  std::vector<arma::uword> perm(std::max<arma::uword>(kv, 1u));

  // --------------------------------------------------------------------------
  // Local loop: estimate SY (variable part) and XCw_cube (for partialling-out),
  // and optionally fill TS and Shat for variable-part hat operator Sv.
  // --------------------------------------------------------------------------
  for (arma::uword z=0; z<nTP; z++) {
    const int*    id = idxT.colptr(z);
    const double* wz = WdT.colptr(z);

    arma::uword nl = 0;
    for (arma::uword j=0; j<NN; j++) {
      const int    k1 = id[j];               // 1-based; NA_INTEGER < 1
      const double w  = wz[j];
      if (k1 >= 1 && arma::uword(k1) <= n && w > 0.0) {
        rows[nl] = arma::uword(k1) - 1u;
        sw[nl]   = std::sqrt(w);
        ++nl;
      }
    }

    if (nl <= kv) continue;

    Xw.set_size(nl, kv);
    XCw.set_size(nl, kc);
    yw.set_size(nl);
    for (arma::uword j=0; j<nl; j++) {
      const double  s  = sw[j];
      const double* xv = XVt.colptr(rows[j]);
      const double* xc = XCt.colptr(rows[j]);
      for (arma::uword c=0; c<kv; c++) Xw.at(j, c)  = xv[c] * s;
      for (arma::uword c=0; c<kc; c++) XCw.at(j, c) = xc[c] * s;
      yw[j] = y[rows[j]] * s;
    }

    // ---- column-pivoted QR of Xw (see gwr_beta_pivotal_qrp_core), Q_r formed
    A = Xw;
    blas_int m_i = blas_int(nl), p_i = blas_int(kv), lda = blas_int(nl);
    blas_int lwork_i = blas_int(work.n_elem), info = 0;
    for (arma::uword c = 0; c < kv; ++c) jpvt[c] = 0;
    arma::lapack::geqp3(&m_i, &p_i, A.memptr(), &lda, jpvt.data(), tau.memptr(),
                        work.memptr(), &lwork_i, &info);
    if (info != 0) continue;
    R = arma::trimatu(A.rows(0, kv - 1));
    const arma::uword r = eff_rank_from_R(R);
    if (r == 0) continue;
    const arma::mat Rr = R.submat(0, 0, r - 1, r - 1);
    for (arma::uword c = 0; c < kv; ++c) perm[c] = arma::uword(jpvt[c]) - 1u;
    blas_int r_i = blas_int(r);
    arma::lapack::orgqr(&m_i, &r_i, &r_i, A.memptr(), &lda, tau.memptr(),
                        work.memptr(), &lwork_i, &info);
    if (info != 0) continue;
    const arma::mat Qr = A.cols(0, r - 1);               // nl x r

    // Local beta for variable part (pivoted basis, then back to the kv columns)
    const arma::vec Qt_y = Qr.t() * yw;
    arma::vec br;
    try {
      RcerrSilencer quiet(sink.rdbuf());
      br = arma::solve(arma::trimatu(Rr), Qt_y, arma::solve_opts::fast);
    } catch (...) {
      continue;
    }
    arma::vec bSY(kv, arma::fill::zeros);
    for (arma::uword k = 0; k < r; ++k) bSY[perm[k]] = br[k];
    SY.row(z) = bSY.t();

    // Local mapping for XC (needed for ZXc = XC - XV * XCwi)
    const arma::mat Qt_Xc = Qr.t() * XCw;
    const arma::mat XCw_r = arma::solve(arma::trimatu(Rr), Qt_Xc, arma::solve_opts::fast);
    for (arma::uword k = 0; k < r; ++k) XCw_cube.slice(z).row(perm[k]) = XCw_r.row(k);

    if (get_se) {
      // Local sigma^2 and SE for variable coefficients (approx local)
      arma::vec yhat_local(nl, arma::fill::zeros);
      for (arma::uword k = 0; k < r; ++k) yhat_local += Xw.col(perm[k]) * br[k];
      const arma::vec res_local = yw - yhat_local;
      const double sigma2_local = arma::dot(res_local, res_local) /
        std::max(1.0, double(nl) - double(r));
      const arma::mat Rinverse = arma::inv(arma::trimatu(Rr));
      const arma::vec sev_loc = arma::sqrt(arma::sum(arma::square(Rinverse), 1) * sigma2_local);
      arma::vec sev_full(kv, arma::fill::zeros);
      for (arma::uword k = 0; k < r; ++k) sev_full[perm[k]] = sev_loc[k];
      SEV.row(z) = sev_full.t();
    }

    // Variable-part hat operator Sv at the focal point: S(f, j) = (sw_j / sw_f) q_f' q_j
    if (get_ts || get_Rk || get_s) {
      const arma::uword focal_id = TP[z] - 1u;
      arma::uword jf = nl;
      for (arma::uword jj = 0; jj < nl; jj++)
        if (rows[jj] == focal_id) { jf = jj; break; }
      if (jf < nl) {
        const arma::vec qf = Qr.row(jf).t();
        arma::vec s = Qr * qf;
        for (arma::uword j = 0; j < nl; j++) s[j] *= sw[j] / sw[jf];
        TS[z] = s[jf];                                   // diag(Sv) at focal = ||q_f||^2
        if (get_s)
          for (arma::uword j = 0; j < nl; j++) Shat(z, rows[j]) = s[j];
      }
    }
  }
  
  // We keep "tS_var" for variable part only
  double tS_var = arma::sum(TS);
  
  // --------------------------------------------------------------------------
  // Fixed part: build ZY and ZXc = (I - Sv) XC at target points
  // --------------------------------------------------------------------------
  arma::vec ZY(nTP, arma::fill::zeros);
  arma::mat ZXc(nTP, kc, arma::fill::zeros);
  
  for (arma::uword z=0; z<nTP; z++) {
    arma::uword i = TP[z] - 1u;

    const arma::mat& XCwi = XCw_cube.slice(z);

    // y* = y - XV * beta_v
    ZY[z] = y[i] - arma::as_scalar(XV.row(i) * SY.row(z).t());
    
    // Z = XC - XV * (mapping of XC through variable local fit)
    ZXc.row(z) = XC.row(i) - XV.row(i) * XCwi;
  }
  
  // beta_c on Z
  arma::mat XtX = ZXc.t() * ZXc;
  arma::vec Xty = ZXc.t() * ZY;
  // the constant coefficients are identified by Z = (I - Sv) XC: a fixed
  // variable that is constant, duplicated or carried by the varying part
  // leaves Z'Z singular, and the former solve() returned huge coefficients
  bool absorbed = false;
  for (arma::uword c = 0; c < kc; ++c) {
    double z2 = 0.0, x2 = 0.0;
    for (arma::uword z = 0; z < nTP; ++z) { z2 += ZXc(z, c) * ZXc(z, c); x2 += XC(TP[z] - 1u, c) * XC(TP[z] - 1u, c); }
    if (!(z2 > 1e-10 * x2)) absorbed = true;   // the varying part reproduces this fixed column
  }
  if (absorbed || !XtX.is_finite() || arma::rcond(arma::symmatu(XtX)) < 1e-12)
    Rcpp::stop("the constant coefficients (fixed_vars) are not identifiable: after removing the varying part, their design is singular (a fixed variable is constant, duplicated, or collinear with the varying coefficients).");
  
  arma::vec beta_c;
  {
    bool ok = arma::solve(beta_c, XtX, Xty, arma::solve_opts::fast);
    if (!ok) {
      // fallback
      beta_c = arma::solve(XtX, Xty);
    }
  }
  
  // Betav = SY - XCwi * beta_c
  arma::mat Betav(nTP, kv, arma::fill::zeros);
  for (arma::uword z=0; z<nTP; z++) {
    Betav.row(z) = SY.row(z) - (XCw_cube.slice(z) * beta_c).t();
  }
  
  // --------------------------------------------------------------------------
  // Compute diag(Hfix) where Hfix = Z (Z'Z)^{-1} Z'  with Z=ZXc
  // and approx diag(S_tot) = diag(Sv) + diag(Hfix)*(1 - diag(Sv))
  // --------------------------------------------------------------------------
  arma::mat XtX_inv;
  {
    bool ok_inv = arma::inv_sympd(XtX_inv, XtX);
    if (!ok_inv) XtX_inv = arma::pinv(XtX);
  }
  
  arma::vec diag_Hfix(nTP, arma::fill::zeros);
  {
    arma::mat Z_A = ZXc * XtX_inv;        // nTP x kc
    diag_Hfix = arma::sum(Z_A % ZXc, 1);  // nTP
  }
  
  // diag(Sv) at focal: prefer Shat if available (guarantees consistency)
  arma::vec diag_Sv = TS;
  if (get_s) {
    for (arma::uword z=0; z<nTP; z++) {
      arma::uword i = TP[z] - 1u;
      diag_Sv[z] = Shat(z, i);
    }
  }
  
  arma::vec diag_S_tot_approx = diag_Sv + diag_Hfix % (1.0 - diag_Sv);
  double tS_tot_approx = arma::sum(diag_S_tot_approx);
  
  // --------------------------------------------------------------------------
  // Global residuals for diagnostics
  // --------------------------------------------------------------------------
  arma::vec resid_full(nTP);
  for (arma::uword z=0; z<nTP; z++) {
    arma::uword i = TP[z] - 1u;
    double fit = arma::as_scalar(XC.row(i) * beta_c)
      + arma::as_scalar(XV.row(i) * Betav.row(z).t());
    resid_full[z] = y[i] - fit;
  }
  
  // Use tS_tot_approx in denominator (recommended for mixed model)
  double denom = double(nTP) - tS_tot_approx - double(kc);
  double sigma2_global = arma::dot(resid_full, resid_full) / std::max(1.0, denom);
  
  arma::vec SE_beta_c;
  if (get_se) {
    arma::vec resid_c = ZY - ZXc * beta_c;
    double sigma2_c = arma::dot(resid_c, resid_c) / std::max(1.0, double(nTP) - double(kc));
    
    // SE(beta_c) = sqrt( sigma2_c * diag((Z'Z)^{-1}) )
    SE_beta_c = arma::sqrt(sigma2_c * XtX_inv.diag());
  }
  
  // --------------------------------------------------------------------------
  // Output
  // --------------------------------------------------------------------------
  Rcpp::List out;
  out["Betav"] = Betav;
  out["Betac"] = beta_c;
  out["SEV"]   = SEV;
  if (get_se) out["se"] = SE_beta_c;
  
  // Backward compatible key:
  out["tS"] = tS_var;
  
  // New diagnostics:
  out["tS_var"]         = tS_var;
  out["tS_tot_approx"]  = tS_tot_approx;
  out["diag_Hfix"]      = diag_Hfix;
  out["diag_Sv"]        = diag_Sv;
  out["diag_S_tot_approx"] = diag_S_tot_approx;
  
  out["sigma2_diag"] = sigma2_global;
  
  if (get_s)  out["Shat"] = Shat;
  if (get_ts) out["TS"]   = TS;
  
  return out;
}
// ---------------------------------------------------------------------------
// Functions exported to R (direct Rcpp interface)
// ---------------------------------------------------------------------------


Rcpp::List gwr_beta_pivotal_qrp_cpp(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericVector& y,
                                    const Rcpp::NumericMatrix& XV,
                                    const Rcpp::IntegerMatrix& indexG,
                                    const Rcpp::NumericMatrix& Wd,
                                    const Rcpp::IntegerVector& TP,
                                    bool get_ts,
                                    bool get_s,
                                    bool get_Rk,
                                    bool get_se) {

  if (Wd.ncol() < indexG.ncol() || Wd.nrow() < TP.size() || indexG.nrow() < TP.size())
    Rcpp::stop("indexG and Wd must have one row per target point and matching columns.");

  arma::mat Xt = arma::mat(const_cast<double*>(X.begin()),
                           X.nrow(), X.ncol(), false).t();
  arma::vec ay(const_cast<double*>(y.begin()), y.size(), false);
  arma::mat aXV(const_cast<double*>(XV.begin()), XV.nrow(), XV.ncol(), false);
  arma::Mat<int> idxT = arma::Mat<int>(const_cast<int*>(indexG.begin()),
                                       indexG.nrow(), indexG.ncol(), false).t();
  arma::mat WdT = arma::mat(const_cast<double*>(Wd.begin()),
                            Wd.nrow(), Wd.ncol(), false).t();
  arma::uvec aTP = Rcpp::as<arma::uvec>(TP);

  return gwr_beta_pivotal_qrp_core(Xt, ay, aXV, idxT, WdT, aTP,
                                   get_ts, get_s, get_Rk, get_se);
}

// -----------------------------------------------------------------------------
// Wrapper: mgwr_beta_pivotal_qrp_mixed_cpp
// -----------------------------------------------------------------------------
Rcpp::List mgwr_beta_pivotal_qrp_mixed_cpp(
    const Rcpp::NumericMatrix& XV,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& XC, // XC is passed as NumericMatrix for attributes
    const Rcpp::IntegerMatrix& indexG,
    const Rcpp::NumericMatrix& Wd,
    const Rcpp::IntegerVector& TP,
    bool get_ts,
    bool get_s,
    bool get_Rk,
    bool get_se
) {
  if (Wd.ncol() < indexG.ncol() || Wd.nrow() < TP.size() || indexG.nrow() < TP.size())
    Rcpp::stop("indexG and Wd must have one row per target point and matching columns.");

  // Views on the R memory, plus transposed copies for contiguous neighbour access
  arma::mat aXV(const_cast<double*>(XV.begin()), XV.nrow(), XV.ncol(), false);
  arma::mat aXC(const_cast<double*>(XC.begin()), XC.nrow(), XC.ncol(), false);
  arma::vec ay(const_cast<double*>(y.begin()), y.size(), false);
  arma::mat XVt = aXV.t();
  arma::mat XCt = aXC.t();
  arma::Mat<int> idxT = arma::Mat<int>(const_cast<int*>(indexG.begin()),
                                       indexG.nrow(), indexG.ncol(), false).t();
  arma::mat WdT = arma::mat(const_cast<double*>(Wd.begin()),
                            Wd.nrow(), Wd.ncol(), false).t();
  arma::uvec aTP = Rcpp::as<arma::uvec>(TP);

  // Call computational kernel
  Rcpp::List out = mgwr_beta_pivotal_qrp_mixed_core_new(
    aXV, XVt, ay, aXC, XCt, idxT, WdT, aTP,
    get_ts, get_s, get_Rk, get_se
  );

  // Handle names for "se" if necessary (as requested in your inline code)
  if (get_se && out.containsElementNamed("se")) {
    Rcpp::NumericVector se_out = Rcpp::wrap(out["se"]); // Wrap the arma::vec

    // Extract names from initial R object XC
    Rcpp::List dn = XC.attr("dimnames");
    if (dn.size() >= 2) {
      Rcpp::CharacterVector xc_names = dn[1];
      se_out.attr("names") = xc_names;
    }
    // Update output list
    out["se"] = se_out;
  }

  return out;
}

// ---------------------------------------------------------------------------
// Kernel weights and row normalisation (native paths of prep_w / normW)
// ---------------------------------------------------------------------------
// kernel_w_cpp reproduces, value for value, the R kernels bisq, gauss, epane
// and triangle and their *_adapt_sorted versions (same formulas, same
// operation order, attributes of d kept). It returns R_NilValue whenever the
// input is outside the domain it reproduces exactly (other kernels, missing
// or negative distances, non-positive bandwidth, column out of range, h of
// unusual length); the caller then falls back to the R kernel.

namespace {

enum KernelCode { K_BISQ, K_GAUSS, K_EPANE, K_TRIANGLE, K_NONE };

// x * x is stored through a volatile so that the compiler cannot fuse it with
// the following subtraction (FMA contraction): R rounds x^2 first.
inline double sq_rounded(double x) {
  volatile double x2 = x * x;
  return x2;
}

template <KernelCode K> inline double kval(double x);
template <> inline double kval<K_BISQ>(double x) {
  if (x < 1.0) { const double t = 1.0 - sq_rounded(x); return (15.0 / 16.0) * (t * t); }
  return 0.0;
}
template <> inline double kval<K_GAUSS>(double x)    { return std::exp(-0.5 * (x * x)); }
template <> inline double kval<K_EPANE>(double x)    { return (x < 1.0) ? (3.0 / 4.0) * (1.0 - sq_rounded(x)) : 0.0; }
template <> inline double kval<K_TRIANGLE>(double x) { return (x <= 1.0) ? (1.0 - x) : 0.0; }

// Column-major fill with raw pointers; one specialised loop per kernel.
template <KernelCode K>
void fill_kernel(const double* dp, double* wp, const double* hr, R_xlen_t nr, R_xlen_t nc) {
  for (R_xlen_t j = 0; j < nc; ++j) {
    const double* dc = dp + j * nr;
    double*       wc = wp + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i) wc[i] = kval<K>(dc[i] / hr[i]);
  }
}

// In-place x / rowSums(x, na.rm = TRUE) with zero sums set to 1: the
// operations of normW()/normw_dense_cpp (row sums over columns in long double).
inline void normalize_rows_inplace(double* xp, R_xlen_t nr, R_xlen_t nc) {
  std::vector<long double> acc(nr, 0.0L);
  for (R_xlen_t j = 0; j < nc; ++j) {
    const double* xc = xp + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i)
      if (!ISNAN(xc[i])) acc[i] += xc[i];
  }
  std::vector<double> rs(nr);
  for (R_xlen_t i = 0; i < nr; ++i) {
    rs[i] = double(acc[i]);
    if (rs[i] == 0.0) rs[i] = 1.0;
  }
  for (R_xlen_t j = 0; j < nc; ++j) {
    double* xc = xp + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i) xc[i] = xc[i] / rs[i];
  }
}

// New double vector with the attributes (dim, dimnames, ...) of `like`,
// without zero-initialisation.
inline Rcpp::NumericVector alloc_like(SEXP like, R_xlen_t len) {
  Rcpp::NumericVector out(Rcpp::no_init(len));
  DUPLICATE_ATTRIB(out, like);
  return out;
}

} // namespace

// n_norm row normalisations (0, 1 or 2) are applied in place, which equals
// normW() applied n_norm times to the kernel matrix, with a single allocation.
SEXP kernel_w_cpp(const Rcpp::NumericMatrix& d, const Rcpp::NumericVector& h,
                  const std::string& kernel, const int n_norm) {
  static const std::string suffix = "_adapt_sorted";
  const bool adaptive = kernel.size() > suffix.size() &&
    kernel.compare(kernel.size() - suffix.size(), suffix.size(), suffix) == 0;
  const std::string base = adaptive ? kernel.substr(0, kernel.size() - suffix.size()) : kernel;

  KernelCode k = K_NONE;
  if (base == "bisq") k = K_BISQ;
  else if (base == "gauss") k = K_GAUSS;
  else if (base == "epane") k = K_EPANE;
  else if (base == "triangle") k = K_TRIANGLE;
  if (k == K_NONE) return R_NilValue;

  const R_xlen_t nr = d.nrow(), nc = d.ncol(), nh = h.size();
  if (nr == 0 || nc == 0 || nh == 0) return R_NilValue;
  const R_xlen_t len = nr * nc;
  const double* dp = d.begin();
  const double* hp = h.begin();

  for (R_xlen_t i = 0; i < len; ++i) {
    const double v = dp[i];
    if (!R_FINITE(v) || v < 0.0) return R_NilValue;
  }

  // Per-row bandwidth
  std::vector<double> hr(nr);
  bool h_const = true;
  for (R_xlen_t i = 1; i < nh; ++i) if (hp[i] != hp[0]) { h_const = false; break; }
  if (!h_const && nh != nr) return R_NilValue;

  if (!adaptive) {
    for (R_xlen_t i = 0; i < nr; ++i) hr[i] = h_const ? hp[0] : hp[i];
  } else {
    // column h (gauss) or h + 2 (compact kernels) of the row-sorted distances
    const double off = (k == K_GAUSS) ? 0.0 : 2.0;
    bool first_row_sorted = true;
    for (R_xlen_t j = 1; j < nc; ++j)
      if (dp[j * nr] < dp[(j - 1) * nr]) { first_row_sorted = false; break; }
    std::vector<double> row(nc);
    for (R_xlen_t i = 0; i < nr; ++i) {
      const double hi = h_const ? hp[0] : hp[i];
      if (!R_FINITE(hi)) return R_NilValue;
      const double cd = hi + off;                  // R truncates double indices
      if (!(cd >= 1.0) || cd >= double(nc) + 1.0) return R_NilValue;
      const R_xlen_t c = R_xlen_t(cd) - 1;
      if (first_row_sorted) {
        hr[i] = dp[i + c * nr];
      } else {
        for (R_xlen_t j = 0; j < nc; ++j) row[j] = dp[i + j * nr];
        std::nth_element(row.begin(), row.begin() + c, row.end());
        hr[i] = row[c];
      }
    }
  }
  for (R_xlen_t i = 0; i < nr; ++i)
    if (!R_FINITE(hr[i]) || !(hr[i] > 0.0)) return R_NilValue;

  Rcpp::NumericVector w = alloc_like(d, len);      // keeps dim and dimnames
  double* wp = w.begin();
  switch (k) {
  case K_BISQ:     fill_kernel<K_BISQ>(dp, wp, hr.data(), nr, nc); break;
  case K_GAUSS:    fill_kernel<K_GAUSS>(dp, wp, hr.data(), nr, nc); break;
  case K_EPANE:    fill_kernel<K_EPANE>(dp, wp, hr.data(), nr, nc); break;
  case K_TRIANGLE: fill_kernel<K_TRIANGLE>(dp, wp, hr.data(), nr, nc); break;
  default: return R_NilValue;
  }
  for (int r = 0; r < n_norm; ++r) normalize_rows_inplace(wp, nr, nc);
  return w;
}

// Dense branch of normW(): x / rowSums(x, na.rm = TRUE), zero sums set to 1.
// Row sums are accumulated column by column in long double, as base rowSums.
SEXP normw_dense_cpp(const Rcpp::NumericMatrix& x) {
  const R_xlen_t nr = x.nrow(), nc = x.ncol();
  Rcpp::NumericVector out = alloc_like(x, nr * nc); // keeps attributes
  std::copy(x.begin(), x.end(), out.begin());
  normalize_rows_inplace(out.begin(), nr, nc);
  return out;
}

// Fused product and row normalisation for GDT weights (alpha = 1):
// normW(Ws * Wt) without materialising the product. Same operations as R:
// double product, row sums over columns in long double skipping NaN, zero
// sums set to 1, division. Attributes are taken from Ws (the caller checks
// that Ws and Wt carry identical attributes).
SEXP wprod_norm_cpp(const Rcpp::NumericMatrix& Ws, const Rcpp::NumericMatrix& Wt) {
  const R_xlen_t nr = Ws.nrow(), nc = Ws.ncol();
  if (Wt.nrow() != nr || Wt.ncol() != nc) Rcpp::stop("Ws and Wt must have the same dimensions.");
  const double* sp = Ws.begin();
  const double* tp = Wt.begin();
  std::vector<long double> acc(nr, 0.0L);
  for (R_xlen_t j = 0; j < nc; ++j) {
    const double* sc = sp + j * nr;
    const double* tc = tp + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i) {
      const double v = sc[i] * tc[i];
      if (!ISNAN(v)) acc[i] += v;
    }
  }
  std::vector<double> rs(nr);
  for (R_xlen_t i = 0; i < nr; ++i) {
    rs[i] = double(acc[i]);
    if (rs[i] == 0.0) rs[i] = 1.0;
  }
  Rcpp::NumericVector out = alloc_like(Ws, nr * nc);
  double* op = out.begin();
  for (R_xlen_t j = 0; j < nc; ++j) {
    const double* sc = sp + j * nr;
    const double* tc = tp + j * nr;
    double*       oc = op + j * nr;
    for (R_xlen_t i = 0; i < nr; ++i) oc[i] = (sc[i] * tc[i]) / rs[i];
  }
  return out;
}
