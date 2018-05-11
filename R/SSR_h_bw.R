#' SSR_h_bw
#' to be documented
#' @usage SSR_h_bw(formula, data,  coord, fixed_vars,kernels, H,
#' h_w,kernel_w='bisq',W1=NULL, Model, control=list(),Penalized)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param kernels to be documented
#' @param H to be documented
#' @param h_w to be documented
#' @param kernel_w to be documented
#' @param W1 to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param Penalized to be documented
#' @keywords internal
#' @return to be documented
SSR_h_bw <-
function (formula, data,  coord, fixed_vars,kernels, H,h_w,kernel_w='bisq',W1=NULL, Model, control=list(),Penalized)
{	#control2=control
	if(is.null(W1) & !(Model %in% c('OLS','GWR','MGWR')))
	  	  control$W<-kernelW_C(coord,h_w,kernel_w,FALSE,Type,0,500,5000,FALSE,0,TRUE) #control$Type --> Type
	SSR_h(formula, data,  coord, fixed_vars,kernels, H, Model, control,Penalized)
}
