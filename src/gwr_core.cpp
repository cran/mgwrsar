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

// ---------------------------------------------------------------------------
// Internal C++ function: gwr_beta_univar_core()
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Internal C++ function: gwr_beta_univar_core()
// ---------------------------------------------------------------------------
Rcpp::List gwr_beta_univar_core(const arma::vec& y,
                                const arma::vec& x,
                                const arma::mat& XV,
                                const arma::umat& indexG,
                                const arma::mat& Wd,
                                const arma::uvec& TP,
                                const bool get_ts,
                                const bool get_s,
                                const bool get_se) {

  const arma::uword nTP = TP.n_elem;
  const arma::uword n   = x.n_elem;

  arma::vec Betav(nTP, arma::fill::zeros);
  arma::vec SEV(nTP, arma::fill::zeros);
  arma::vec TS(nTP, arma::fill::zeros);
  arma::mat Shat;
  if (get_s) Shat = arma::zeros(nTP, n);

  for (arma::uword z = 0; z < nTP; ++z) {
    arma::uvec idx = indexG.row(z).t() - 1u;
    arma::vec wi   = Wd.row(z).t();

    std::vector<arma::uword> keep;
    std::vector<double> wi_vec;
    keep.reserve(idx.n_elem);
    wi_vec.reserve(idx.n_elem);
    for (arma::uword j = 0; j < idx.n_elem; ++j)
      if (wi[j] > 1e-12 && idx[j] < n) {
        keep.push_back(idx[j]);
        wi_vec.push_back(wi[j]);
      }

      if (keep.empty()) continue;

      arma::uvec rows = arma::conv_to<arma::uvec>::from(keep);
      arma::vec x_loc = x.elem(rows);
      arma::vec y_loc = y.elem(rows);
      arma::vec w_loc = arma::conv_to<arma::vec>::from(wi_vec);

      double Sxx = arma::dot(w_loc % x_loc, x_loc);
      double Sxy = arma::dot(w_loc % x_loc, y_loc);
      double beta = (Sxx == 0.0) ? 0.0 : Sxy / Sxx;
      Betav[z] = beta;

      // --- Compute SEV (must be INSIDE the z loop) ---
      if (get_se) {
        if (Sxx > 0.0) {
          arma::vec res_loc = (y_loc - beta * x_loc);
          double rss_loc = arma::dot(w_loc % res_loc, res_loc);
          double sigma2_loc = rss_loc / std::max(1.0, double(rows.n_elem) - 1.0);
          SEV[z] = std::sqrt(sigma2_loc / Sxx);
        }
      }

      if (get_ts || get_s) {
        double x0 = XV(TP[z] - 1u, 0);
        arma::vec s_local(rows.n_elem, arma::fill::zeros);
        if (Sxx != 0.0)
          s_local = w_loc % (x_loc * (x0 / Sxx));

        if (get_ts && s_local.n_elem > 0) {
          arma::uword focal_id = TP[z] - 1u;
          for (arma::uword jj = 0; jj < rows.n_elem; jj++) {
            if (rows[jj] == focal_id) {
              TS[z] = s_local[jj];
              break;
            }
          }
        }

        if (get_s)
          for (arma::uword j = 0; j < rows.n_elem; ++j)
            Shat(z, rows[j]) = s_local[j];
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

  arma::vec ay(const_cast<double*>(y.begin()), y.size(), false);
  arma::vec ax(const_cast<double*>(x.begin()), x.size(), false);
  arma::mat aXV(const_cast<double*>(XV.begin()), XV.nrow(), XV.ncol(), false);
  arma::umat aIndexG = Rcpp::as<arma::umat>(indexG);
  arma::mat aWd(const_cast<double*>(Wd.begin()), Wd.nrow(), Wd.ncol(), false);
  arma::uvec aTP = Rcpp::as<arma::uvec>(TP);

  // Call with all 9 arguments
  return gwr_beta_univar_core(ay, ax, aXV, aIndexG, aWd, aTP, get_ts, get_s, get_se);
}
// ---------------------------------------------------------------------------
// Main function: gwr_beta_pivotal_qrp_core (multivariate pivoted QR)
// ---------------------------------------------------------------------------

Rcpp::List gwr_beta_pivotal_qrp_core(
    const arma::mat& X,
    const arma::vec& y,
    const arma::mat& XV,
    const arma::umat& indexG,
    const arma::mat& Wd,
    const arma::uvec& TP,
    const bool get_ts,
    const bool get_s,
    const bool get_Rk,
    const bool get_se
) {
  using namespace arma;

  const uword nTP = TP.n_elem;
  const uword n   = X.n_rows;
  const uword p   = X.n_cols;

  mat Betav(nTP, p, fill::zeros);

  vec TS(nTP, fill::zeros);
  mat Shat; if (get_s)  Shat = mat(nTP, n, fill::zeros);
  cube Rk;  if (get_Rk) Rk   = cube(nTP, n, p, fill::zeros);

  mat SEV;  if (get_se) SEV  = mat(nTP, p, fill::zeros);

  // ---- loop over focal points
  for (uword z = 0; z < nTP; ++z) {

    ivec idx_i = conv_to<ivec>::from(indexG.row(z).t()) - 1;
    rowvec wi_full = Wd.row(z);

    std::vector<uword> rows_vec;
    std::vector<double> wi_vec;
    rows_vec.reserve(idx_i.n_elem);
    wi_vec.reserve(idx_i.n_elem);

    // keep only valid neighbors (and keep their weights aligned)
    for (uword j = 0; j < idx_i.n_elem; ++j) {
      int id = idx_i[j];
      double w = wi_full[j];
      if (id >= 0 && id < (int)n && w > 0) {
        rows_vec.push_back((uword)id);
        wi_vec.push_back(w);
      }
    }

    const uword n_loc = rows_vec.size();
    if (n_loc <= p) continue;

    uvec rows = conv_to<uvec>::from(rows_vec);
    vec  wi   = conv_to<vec>::from(wi_vec);

    vec sw = sqrt(wi);

    mat X_loc = X.rows(rows);
    vec y_loc = y.elem(rows);
    if (X_loc.n_rows == 0) continue;

    mat Xw = X_loc.each_col() % sw;
    vec yw = y_loc % sw;

    mat Q, R;
    try { qr_econ(Q, R, Xw); } catch (...) { continue; }

    uword r = eff_rank_from_R(R);
    if (r == 0u) continue;

    mat Rr = R.submat(0, 0, r - 1, r - 1);
    vec Qt_y = Q.t() * yw;

    vec beta_r;
    try {
      std::stringstream buffer;
      std::streambuf* old_buf = Rcpp::Rcerr.rdbuf(buffer.rdbuf());
      beta_r = solve(trimatu(Rr), Qt_y.head(r), solve_opts::fast);
      Rcpp::Rcerr.rdbuf(old_buf);
    } catch (...) {
      continue;
    }

    // expand to full p with zeros on non-estimable columns
    vec beta(p, fill::zeros);
    beta.head(r) = beta_r;
    Betav.row(z) = beta.t();

    // ---- local standard errors (only if requested)
    if (get_se) {
      vec sev_full(p, fill::zeros);

      try {
        // local fitted values and weighted residuals
        vec yhat_w = Xw.cols(0, r - 1) * beta_r; // because columns beyond r are zero
        vec resid_w = yw - yhat_w;

        double rss_w = dot(resid_w, resid_w);
        double denom = std::max(1.0, double(n_loc) - double(r));
        double sigma2 = rss_w / denom;

        mat Rinverse = inv(trimatu(Rr));        // r x r
        mat XtXinv_r = Rinverse * Rinverse.t(); // r x r

        vec sev_loc = sqrt(sigma2 * XtXinv_r.diag());

        for (uword j = 0; j < r; ++j) {
          double v = sev_loc[j];
          sev_full[j] = (std::isfinite(v) ? v : 0.0);
        }
      } catch (...) {
        // keep zeros (legacy behavior)
      }

      SEV.row(z) = sev_full.t();
    }

    // ---- TS / Shat / Rk (hat-related outputs)
    if (get_ts || get_s || get_Rk) {

      const rowvec x0 = XV.row(TP[z] - 1u);
      mat QrT = Q.cols(0, r - 1).t();
      mat Br;

      try {
        std::stringstream buffer;
        std::streambuf* old_buf = Rcpp::Rcerr.rdbuf(buffer.rdbuf());
        Br = solve(trimatu(Rr), QrT, solve_opts::fast);
        Rcpp::Rcerr.rdbuf(old_buf);
      } catch (...) {
        continue;
      }

      mat B(p, n_loc, fill::zeros);
      B.rows(0, r - 1) = Br;

      vec u = (x0 * B).t();
      vec s_local = sw % u;


      if (get_ts && s_local.n_elem > 0) {
        uword focal_id = TP[z] - 1u;
        TS[z] = 0.0;
        for (uword jj = 0; jj < n_loc; jj++) {
          if (rows[jj] == focal_id) {
            TS[z] = s_local[jj];
            break;
          }
        }
      }

      if (get_s)
        for (uword j = 0; j < n_loc; ++j)
          Shat(z, rows[j]) = s_local[j];

      if (get_Rk) {
        const uword focal = TP[z] - 1u;
        for (uword nx = 0; nx < p; ++nx) {
          double x0nx = x0[nx];
          if (!std::isnan(x0nx)) {
            for (uword j = 0; j < n_loc; ++j)
              Rk(focal, rows[j], nx) = x0nx * B(nx, j) * sw[j];
          }
        }
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
// ============================================================================

Rcpp::List mgwr_beta_pivotal_qrp_mixed_core_new(
    const arma::mat& XV,
    const arma::vec& y,
    const arma::mat& XC,
    const arma::umat& indexG,
    const arma::mat& Wd,
    const arma::uvec& TP,
    const bool get_ts,
    const bool get_s,
    const bool get_Rk,
    const bool get_se
) {
  
  arma::uword nTP = TP.n_elem,
    n   = XV.n_rows,
    kv  = XV.n_cols,
    kc  = XC.n_cols;
  
  arma::mat SY(nTP, kv, arma::fill::zeros);
  arma::mat XCw_flat(nTP * kv, kc, arma::fill::zeros);
  
  // TS stores diag(Sv) at focal point (i=TP[z]-1) for each z
  arma::vec TS(nTP, arma::fill::zeros);
  
  arma::mat SEV(nTP, kv, arma::fill::zeros);
  
  arma::mat Shat;
  if (get_s) Shat = arma::zeros(nTP, n);
  
  arma::cube Rk;
  if (get_Rk) Rk = arma::cube(n, n, kv, arma::fill::zeros); // (unused here, kept for API)
  
  // --------------------------------------------------------------------------
  // Local loop: estimate SY (variable part) and XCw_flat (for partialling-out),
  // and optionally fill TS and Shat for variable-part hat operator Sv.
  // --------------------------------------------------------------------------
  for (arma::uword z=0; z<nTP; z++) {
    arma::ivec idx = arma::conv_to<arma::ivec>::from(indexG.row(z).t()) - 1;
    arma::rowvec wrow = Wd.row(z);
    
    std::vector<arma::uword> rv;
    std::vector<double> wv;
    rv.reserve(idx.n_elem);
    wv.reserve(idx.n_elem);
    
    for (arma::uword j=0; j<idx.n_elem; j++) {
      int id = idx[j];
      double w = wrow[j];
      if (id >= 0 && id < (int)n && w > 0.0) {
        rv.push_back((arma::uword)id);
        wv.push_back(w);
      }
    }
    
    arma::uword nl = rv.size();
    if (nl <= kv) continue;
    
    arma::uvec rows = arma::conv_to<arma::uvec>::from(rv);
    arma::vec wi    = arma::conv_to<arma::vec>::from(wv);
    arma::vec sw    = arma::sqrt(wi);
    
    arma::mat XVloc = XV.rows(rows);
    arma::mat XCloc = XC.rows(rows);
    
    arma::mat Xw  = XVloc.each_col() % sw;
    arma::vec yw  = y.elem(rows) % sw;
    arma::mat XCw = XCloc.each_col() % sw;
    
    arma::mat Q, R;
    arma::qr_econ(Q, R, Xw);
    
    arma::uword r = eff_rank_from_R(R);
    if (r == 0) continue;
    
    arma::mat Rr = R.submat(0, 0, r-1, r-1);
    arma::vec Qt_y = Q.t() * yw;
    
    // Local beta for variable part (in QR basis, then padded to kv)
    arma::vec br = arma::solve(arma::trimatu(Rr), Qt_y.head(r), arma::solve_opts::fast);
    arma::vec bSY(kv, arma::fill::zeros);
    bSY.head(r) = br;
    SY.row(z) = bSY.t();
    
    // Local mapping for XC (needed for ZXc = XC - XV * XCwi)
    arma::mat Qt_Xc = Q.t() * XCw;
    arma::mat XCw_r = arma::solve(arma::trimatu(Rr), Qt_Xc.rows(0, r-1), arma::solve_opts::fast);
    arma::mat XCwi(kv, kc, arma::fill::zeros);
    XCwi.rows(0, r-1) = XCw_r;
    
    arma::uword base = z * kv;
    for (arma::uword rr=0; rr<kv; rr++) {
      XCw_flat.row(base + rr) = XCwi.row(rr);
    }
    
    if (get_se) {
      // Local sigma^2 and SE for variable coefficients (approx local)
      arma::vec yhat_local = Xw.cols(0, r-1) * br;
      arma::vec res_local  = yw - yhat_local;
      double sigma2_local  = arma::dot(res_local, res_local) /
        std::max(1.0, double(nl) - double(r));
      
      arma::mat Rinverse = arma::inv(arma::trimatu(Rr));
      arma::vec sev_loc  = arma::sqrt( arma::diagvec(Rinverse * Rinverse.t()) * sigma2_local );
      
      arma::vec sev_full(kv, arma::fill::zeros);
      sev_full.head(r) = sev_loc.head(r);
      SEV.row(z) = sev_full.t();
    }
    
    // Variable-part hat operator Sv at focal point
    if (get_ts || get_s || get_Rk) {
      
      arma::rowvec x0 = XV.row(TP[z] - 1u);
      
      arma::mat Q_sub = Q.cols(0, r-1);
      arma::mat Br    = arma::solve(arma::trimatu(Rr), Q_sub.t(), arma::solve_opts::fast);
      
      arma::mat B(kv, nl, arma::fill::zeros);
      B.rows(0, r-1) = Br;
      
      arma::vec u = (x0 * B).t();
      arma::vec s = sw % u; // hat row weights at location z, restricted to "rows"
      
      if ((get_ts || get_s || get_Rk) && s.n_elem > 0) {
        arma::uword focal_id = TP[z] - 1u;
        for (arma::uword jj = 0; jj < nl; jj++) {
          if (rows[jj] == focal_id) {
            TS[z] = s[jj]; // diag(Sv) at focal
            break;
          }
        }
      }
      
      if (get_s) {
        for (arma::uword j=0; j<nl; j++) {
          Shat(z, rows[j]) = s[j];
        }
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
    
    arma::mat XCwi(kv, kc);
    for (arma::uword rr=0; rr<kv; rr++) {
      XCwi.row(rr) = XCw_flat.row(z * kv + rr);
    }
    
    // y* = y - XV * beta_v
    ZY[z] = y[i] - arma::as_scalar(XV.row(i) * SY.row(z).t());
    
    // Z = XC - XV * (mapping of XC through variable local fit)
    ZXc.row(z) = XC.row(i) - XV.row(i) * XCwi;
  }
  
  // beta_c on Z
  arma::mat XtX = ZXc.t() * ZXc;
  arma::vec Xty = ZXc.t() * ZY;
  
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
    arma::mat XCwi(kv, kc);
    for (arma::uword rr=0; rr<kv; rr++) {
      XCwi.row(rr) = XCw_flat.row(z * kv + rr);
    }
    Betav.row(z) = SY.row(z) - (XCwi * beta_c).t();
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

  arma::mat aX(const_cast<double*>(X.begin()), X.nrow(), X.ncol(), false);
  arma::vec ay(const_cast<double*>(y.begin()), y.size(), false);
  arma::mat aXV(const_cast<double*>(XV.begin()), XV.nrow(), XV.ncol(), false);
  arma::umat aIndexG = Rcpp::as<arma::umat>(indexG);
  arma::mat aWd(const_cast<double*>(Wd.begin()), Wd.nrow(), Wd.ncol(), false);
  arma::uvec aTP = Rcpp::as<arma::uvec>(TP);

  return gwr_beta_pivotal_qrp_core(aX, ay, aXV, aIndexG, aWd, aTP,
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
  // Clean conversion to Armadillo without unnecessary copy if possible
  arma::mat aXV = Rcpp::as<arma::mat>(XV);
  arma::vec ay  = Rcpp::as<arma::vec>(y);
  arma::mat aXC = Rcpp::as<arma::mat>(XC);
  arma::umat aIndexG = Rcpp::as<arma::umat>(indexG);
  arma::mat aWd = Rcpp::as<arma::mat>(Wd);
  arma::uvec aTP = Rcpp::as<arma::uvec>(TP);

  // Call computational kernel
  Rcpp::List out = mgwr_beta_pivotal_qrp_mixed_core_new(
    aXV, ay, aXC, aIndexG, aWd, aTP,
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
