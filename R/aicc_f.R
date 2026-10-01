#' Corrected Akaike Information Criterion (AICc) for GWR
#'
#' Computes the bias-corrected AIC for models with a smoothing matrix S.
#' Formula: \code{n*log(RSS/n) + n*log(2*pi) + n*(n + pena*tS) / (n - 1 - pena*tS)}.
#'
#' @usage aicc_f(e, ts, n, pena = 1)
#' @param e Numeric vector of residuals.
#' @param ts Trace of the smoothing matrix S (effective number of parameters).
#' @param n Sample size (overridden internally by \code{length(e)}).
#' @param pena Penalty multiplier for the trace, default 1.
#' @noRd
#' @return A scalar AICc value.
aicc_f <-
  function(e,ts,n,pena=1){
    n=length(e)
    # a fit that uses at least n - 1 degrees of freedom is saturated: the
    # penalty term would change sign and reward it instead of rejecting it
    if (n - 1 - pena*ts <= 0) return(Inf)
    #n*log(sum(e^2)/n)+n*log(2*pi)+n*(n+ts)/(n-2-ts)
    n*log(sum(e^2)/n)+n*log(2*pi)+n*(n+pena*ts)/(n-1-pena*ts)

  }
