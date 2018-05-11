#' SSR_h
#' to be documented
#' @usage SSR_h(formula, data,  coord, fixed_vars,kernels, H, Model, control=list(),Penalized)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param kernels to be documented
#' @param H to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param Penalized to be documented
#' @keywords internal
#' @return to be documented
SSR_h <-
function (formula, data,  coord, fixed_vars,kernels, H, Model, control=list(),Penalized)
{
    model <- MGWRSAR(formula, data, coord, fixed_vars, kernels, H, Model,control)
        if(Penalized) corpen=ifelse(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv','MGWRSAR_1_kc_0'),sum(abs(model$Betav[,ncol(model$Betav)])>1),0) else corpen=0
    SSR(model)*(1+corpen)
}
