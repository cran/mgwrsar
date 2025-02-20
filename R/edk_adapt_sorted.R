#' edk_adapt_sorted
#' to be documented
#' @usage gauss_adapt_sorted(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of number of neighbours as bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
edk_adapt_sorted<-function (d, h) {
  if(is.unsorted(d[1,])) {
    dd=t(apply(d,1,sort))
  } else dd=d
  if(any(h!=h[1])) {
    h <- dd[cbind(1:nrow(dd),h)]
  } else h=dd[,h[1]]
  exp(-(d/h))
}
