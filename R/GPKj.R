#' GPKj
#' to be documented
#' @usage GPKj(j, Y,X,S, H, kernels, type = "GD", minv = 0,
#' maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay)
#' @param j  to be documented
#' @param Y  to be documented
#' @param X  to be documented
#' @param S  to be documented
#' @param H  to be documented
#' @param kernels  to be documented
#' @param type  to be documented
#' @param minv  to be documented
#' @param maxknn  to be documented
#' @param NmaxDist  to be documented
#' @param isgcv  to be documented
#' @param TIME  to be documented
#' @param decay  to be documented
#' @keywords internal
#' @return to be documented
GPKj <-
function(j, Y,X,S, H, kernels, type = "GD", minv = 0, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay){
Wd<-do.call(type,args=list(j, Y,X,S, H, kernels, minv = minv, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay) )
Wd <- Wd/sum(Wd)
Wd
}
