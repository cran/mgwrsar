#' cv_h
#' to be documented
#' @usage cv_h(H,Hp,kernel_w=NULL,search_adaptive=TRUE,formula, data,
#' coords, fixed_vars , kernels, Model, control)
#' @param H  A vector of bandwidths
#' @param Hp A vector of bandwidths with NULL for the bandwidth under optimization
#' @param kernel_w The type of kernel for computing W, default NULL
#' @param search_adaptive  A vector of boolean to choose adaptive version for each kernel.
#' @param formula  a formula.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coords default NULL, a dataframe or a matrix with coordinates, not
#' required if data is a spatial dataframe.
#' @param fixed_vars a vector with the names of spatiallay constant coefficient for
#' mixed model. All other variables present in formula are supposed to be spatially
#' varying. If empty or NULL (default), all variables in formula are supposed to be
#' spatially varying.
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"), gaussian
#' ("gauss")).
#' @param Model character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' explanation.
#' @param control list of extra control arguments for MGWRSAR wrapper - see Details below
#' @noRd
cv_h<-function(H,Hp,kernel_w=NULL,search_adaptive=TRUE,formula, data, coords, fixed_vars , kernels, Model, control){
  ### search
  if(!is.null(kernel_w)){
    h_w=H
    if(search_adaptive) NNN=h_w+2 else NNN=500
    control$W=kernel_matW(H=h_w,kernels=kernel_w,coord_i=coords,NN=NNN,Type='GD',adaptive=search_adaptive,diagnull=TRUE,rowNorm=TRUE)
  } else if(length(Hp)>1) Hp[which(is.na(Hp))]<-H else Hp=H
  ###
  model<-MGWRSAR(formula=formula, data = data,coords=coords, fixed_vars=fixed_vars,Model=Model,kernels=kernels,H=Hp,control=control);
  if(control$tp_rmse==1) model$RMSEn else if(control$tp_rmse==2)  sqrt(mean(model$residuals[model$TP]^2)) else if(control$tp_rmse==3) sqrt(mean(model$residuals[-model$TP]^2))
  #sqrt(mean(model$residuals[model$TP]^2))
}
