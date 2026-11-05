#' Saturation filter of the TDS candidate fits
#'
#' A candidate univariate fit whose residual variance falls below \code{ratio}
#' times the residual variance of the current complete fit (the fit before the
#' component update) is a near-interpolating local fit: it is excluded from
#' the arg-min of the AICc over the candidates. The rule is a ratio of
#' variances, so it does not depend on the unit of the response nor on
#' whether the projection matrices are computed; \code{ratio = 0} disables it.
#' Until 1.4.3 the rule compared the AICc of a candidate with half of a
#' reference AICc, which depended on the unit of the response and on the
#' inference option.
#'
#' @param res Matrix of candidate results (one row per candidate, columns
#'   \code{AICc} and \code{e0} among others), as built in \code{update_opt()}
#'   and \code{update_opt_st()}.
#' @param ref_var Residual variance of the current complete fit.
#' @param ratio Threshold on the ratio of residual variances.
#' @return A list: \code{valid} (indices of the admissible candidates),
#'   \code{excluded} (indices of the saturated ones), \code{AICc} and
#'   \code{var_cand} (per candidate).
#' @noRd
.tds_saturation_filter <- function(res, ref_var, ratio) {
  AICc_vals <- as.numeric(res[, "AICc"])
  var_cand <- vapply(seq_len(nrow(res)), function(j) {
    e <- unlist(res[j, "e0"])
    if (!length(e)) return(NA_real_)
    mean(e^2, na.rm = TRUE)
  }, numeric(1))
  ok <- !is.na(AICc_vals)
  saturated <- ok & is.finite(var_cand) & isTRUE(ratio > 0) &
    is.finite(ref_var) & var_cand < ratio * ref_var
  list(valid = which(ok & !saturated), excluded = which(saturated),
       AICc = AICc_vals, var_cand = var_cand)
}
