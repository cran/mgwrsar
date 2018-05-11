#' aic_model
#' to be documented
#' @usage aic_model(model)
#' @param model  to be documented
#' @keywords internal
#' @return to be documented
aic_model <-
function(model){
n=length(model$Y)
log(model$SSR) + 2*(ifelse(is.null(model$tS),model$edf,model$tS))/n
}
