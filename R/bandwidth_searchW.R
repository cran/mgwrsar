#' bandwidth_searchW
#' to be documented
#' @usage bandwidth_searchW(formula,data,coord,fixed_vars,Model,control,kernel,kernel_w,e_search)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param kernel to be documented
#' @param kernel_w to be documented
#' @param e_search to be documented
#' @keywords internal
#' @return to be documented
bandwidth_searchW <-
function(formula,data,coord,fixed_vars,Model,control,kernel,kernel_w,e_search) {
cat('Bandwidth h_w for W: Searching optimal bandwidth for kernel',kernel_w,' with Model',Model,'  ... \n')
if(kernel_w %in% c('bisq','gauss','bin','epanechnikov')) {
if(e_search$lower_cW<e_search$upper_c & e_search$upper_cw>e_search$lower_cW){
opt_bandwidth=try(incr_search_band(
SSR_h_bw,formula=formula, data=data,  coord=coord, fixed_vars=fixed_vars,kernels=kernel,H=e_search$h,  kernel_w = kernel_w, W1 =NULL, Model=Model, control=control, Penalized=e_search$Penalized,discrete=FALSE, lower=e_search$lower_cW,n=e_search$n),silent=TRUE) } else {opt_bandwidth='NOSUPPORT';}
} else  if(kernel_w %in% c('bisq_knn','gauss_adapt','gauss_knn','knn','epanechnikov_knn')){
if(e_search$lower_dW<0.65*e_search$n){
opt_bandwidth=try(incr_search_band(
SSR_h_bw,formula=formula, data=data,  coord=coord, fixed_vars=fixed_vars,kernels=kernel,H=e_search$h,  kernel_w = kernel_w, W1 =NULL, Model=Model, control=control, Penalized=e_search$Penalized,discrete=TRUE, lower=e_search$lower_dW,n=e_search$n),silent=TRUE)
} else {opt_bandwidth='NOSUPPORT';}
}
if(class(opt_bandwidth)=='try-error') opt_bandwidth='NOSUPPORT'
opt_bandwidth
}
