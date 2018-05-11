#' model_conf
#' to be documented
#' @usage model_conf(formula,data,coord,fixed_vars,Model,control,config_model,search_W)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param config_model to be documented
#' @param search_W to be documented
#' @keywords internal
#' @return to be documented
model_conf <-
function(formula,data,coord,fixed_vars,Model,control,config_model,search_W){
if(search_W){
control$kernel_w<-config_model[4]
control$h_w<-as.numeric(config_model[5])
control$W<-kernelW_C(coord,control$h_w,control$kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE) #control$Type --> Type
}
control$isgcv=FALSE
MGWRSAR(formula=formula,data=data,coord=coord,fixed_vars=fixed_vars,kernels=unlist(config_model[2]),H=as.numeric(config_model[3]),Model=Model,control=control)
}
