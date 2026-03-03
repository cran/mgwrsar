#' Proportional AIC for Model Comparison
#'
#' Computes a proportional AIC (divided by n) for relative model comparison.
#' Formula: \code{log(SSR) + 2*tS/n}. Note: this is NOT on the same scale as
#' the \code{@@AIC} slot values computed in \code{format_and_diagno}.
#'
#' @usage aic_model(model)
#' @param model An \code{mgwrsar} S4 model object with slots SSR, tS, edf, Y.
#' @noRd
#' @return A scalar proportional AIC value.
aic_model <-
function(model){
  n=length(model@Y)
log(model@SSR) + 2*(ifelse(is.null(model@tS),model@edf,model@tS))/n
}
