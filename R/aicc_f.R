#' aicc_f
#' to be documented
#' @usage aicc_f(e,ts)
#' @param e  residuals
#' @param ts  trace(S)
#' @noRd
#' @return to be documented
aicc_f <-
  function(e,ts,n,pena=1){
    n=length(e)
    #n*log(sum(e^2)/n)+n*log(2*pi)+n*(n+ts)/(n-2-ts)
    n*log(sum(e^2)/n)+n*log(2*pi)+n*(n+pena*ts)/(n-1-pena*ts)

  }
