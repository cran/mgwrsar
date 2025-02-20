#' fill_DGPTAB
#' to be documented
#' @usage fill_DGPTAB(formula,data,coords,fixed_vars,Model,control,opt,
#' search_W,config_model)
#' @param formula  to be documented
#' @param data  to be documented
#' @param coords  to be documented
#' @param fixed_vars  to be documented
#' @param Model  to be documented
#' @param control  to be documented
#' @param opt  to be documented
#' @param search_W  to be documented
#' @param config_model  to be documented
#' @noRd
#' @return to be documented
fill_DGPTAB <-function(formula,data,coords,fixed_vars,Model,control,opt,search_W,config_model){
  est_model=model_conf(formula,data,coords,fixed_vars,Model,control,config_model,search_W)
  if(!(Model %in% c('OLS','SAR'))){
    list(config_model=config_model,CV=opt$objective,SSR=est_model@SSR,model=est_model,dists=control$dists,indexG=control$indexG)
  } else {list(config_model=config_model,CV=0,SSR=est_model@SSR,model=est_model)}
}
