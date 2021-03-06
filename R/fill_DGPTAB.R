#' fill_DGPTAB
#' to be documented
#' @usage fill_DGPTAB(formula,data,coord,fixed_vars,Model,control,opt,search_W,config_model)
#' @param formula  to be documented
#' @param data  to be documented
#' @param coord  to be documented
#' @param fixed_vars  to be documented
#' @param Model  to be documented
#' @param control  to be documented
#' @param opt  to be documented
#' @param search_W  to be documented
#' @param config_model  to be documented
#' @keywords internal
#' @return to be documented
fill_DGPTAB <-
function(formula,data,coord,fixed_vars,Model,control,opt,search_W,config_model){
    est_model=model_conf(formula,data,coord,fixed_vars,Model,control,config_model,search_W)
    if(!(Model %in% c('OLS','SAR'))){
        list(config_model=config_model,CV=opt$objective,SSR=est_model$SSR,model=est_model)
    } else {list(config_model=config_model,CV=0,SSR=est_model$SSR,model=est_model)} ##(CV(est_model))$CV
}
