#' aic_model
#' This funciton estimate a MGWRSAR model and return the AICc.
#' @usage AICc(H,formula, data ,coords, fixed_vars,kernels, Model,control)
#' @param H vector containing the bandwidth parameters for the kernel functions.
#' @param formula  a formula.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coords default NULL, a dataframe or a matrix with coordinates, not
#' required if data is a spatial dataframe.
#' @param fixed_vars a vector with the names of spatiallay constant coefficient for
#' mixed model. All other variables present in formula are supposed to be spatially
#' varying. If empty or NULL (default), all variables in formula are supposed to be
#' spatially varying.
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"),
#' gaussian ("gauss")) .
#' @param Model character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0".
#' @param control list of extra control arguments for MGWRSAR wrapper - see MGWRSAR
#' @noRd
#' @return AICc value.
AICc<-function(H,formula, data ,coords, fixed_vars,kernels, Model,control){
  control$get_ts=TRUE
  control$NN=H+2
  model<-MGWRSAR(formula=formula , data = data,coords=coords, fixed_vars=fixed_vars,kernels=kernels,H=H, Model = Model,control=control)
  n <- length(model$Y)
  AICc=n*log(model$SSR/n)+n*log(2*pi)+n*(n+model$tS)/(n-2-model$tS)
  #cat('H ',H,' AICc ',AICc,'\n')
  AICc
}
