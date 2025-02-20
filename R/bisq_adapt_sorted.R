#' bisq_adapt_sorted
#' to be documented
#' @usage bisq_adapt_sorted(d,h)
#' @param d a vector of distances
#' @param h a scalar or vector of number of neighbours as bandiwdths (lenght(h)=1 or length(h)=length(d))
#' @noRd
#' @return a vector of weights.
bisq_adapt_sorted <- function(d,h) {     #adaptive bisquare kernel
  if(is.unsorted(d[1,])) {
    dd=t(apply(d,1,sort))
  } else dd=d
  if(any(h!=h[1])) {h <- dd[cbind(1:nrow(dd),h+2)]} else h=dd[,h[1]+2]
  x<-d/h
  x[abs(x)< 1]<- (15/16)*((1-x[x< 1]^2)^2)
  x[abs(x)>= 1]<- 0
  x
}
