#' rectangle
#' to be documented
#' @usage rectangle(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
rectangle<-function (d,h) {
  if(is(d,'matrix')) matrix(as.numeric(d<h),ncol=ncol(d)) else as.numeric(d<h)
}
