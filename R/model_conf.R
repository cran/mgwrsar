#' model_conf
#' to be documented
#' @usage model_conf(formula,data,coords,fixed_vars,Model,control,config_model,search_W)
#' @param formula to be documented
#' @param data to be documented
#' @param coords to be documented
#' @param fixed_vars to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param config_model to be documented
#' @param search_W to be documented
#' @noRd
#' @return to be documented
model_conf <-function(formula,data,coords,fixed_vars,Model,control,config_model,search_W){
  if(search_W){
    control$kernel_w<-config_model$kernel_w
    control$h_w<-as.numeric(config_model$h_w)
  }
  control$isgcv=FALSE
  kernels=unlist(config_model$kernels)
  H=as.numeric(config_model$H)
  MGWRSAR(formula=formula,data=data,coords=coords,fixed_vars=fixed_vars,kernels=kernels,H=H,Model=Model,control=control)
}
