#' bandwidth_search
#' to be documented
#' @usage bandwidth_search(formula,data,coord,fixed_vars,Model,control,kernel,e_search)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param kernel to be documented
#' @param e_search to be documented
#' @keywords internal
#' @return to be documented
bandwidth_search <-
function(formula,data,coord,fixed_vars,Model,control,kernel,e_search){
cat('Bandwidth h: Searching optimal bandwidth for kernel',kernel[1],' with Model',Model,'  ... \n')
if(kernel[1] %in% c('bisq','gauss','bin','epanechnikov','bisq_ext','gauss_ext')) {
if(e_search$lower_cw<e_search$upper_cw & e_search$upper_cw>e_search$lower_cw){
opt_bandwidth=try(optimize(
SSR_h_bw,formula=formula, data=data,  coord=coord, fixed_vars=fixed_vars,kernels=kernel,h_w=e_search$h_w,  kernel_w = e_search$kernels_w[1], W1 =control$W, Model=Model, control=control, Penalized=e_search$Penalized, lower=e_search$lower_cw,upper=e_search$upper_cw),silent=TRUE) } else {opt_bandwidth='NOSUPPORT';}
} else  if(kernel[1] %in% c('bisq_knn','gauss_adapt','gauss_knn','knn','epanechnikov_knn')){
if(e_search$lower_dw<0.65*e_search$n){
opt_bandwidth=try(incr_search_band(
SSR_h_bw,formula=formula, data=data,  coord=coord, fixed_vars=fixed_vars,kernels=kernel,h_w=e_search$h_w,  kernel_w = e_search$kernels_w[1], W1 =control$W, Model=Model, control=control, Penalized=e_search$Penalized,discrete=TRUE, lower=e_search$lower_dw,n=e_search$n),silent=TRUE)
} else {opt_bandwidth='NOSUPPORT';}
}
if(class(opt_bandwidth)=='try-error') opt_bandwidth='NOSUPPORT'
#cat('\n')
opt_bandwidth
}
