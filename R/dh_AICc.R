#' dh_AICc
#' to be documented
#' @usage dh_AICc(e0,e1,ds,wi,i)
#' @param e0 a vector of residus with upper bandwidth
#' @param e1 a vector of residus with lower bandwidth
#' @param ds0 diagonal of Hat Matrix S
#' @param ds1 diagonal of Hat Matrix S
#' @param wi a vector of weigths for i with upper bandwidth
#' @noRd
#' @return a vector of weights.
dh_AICc<-function (e0,e1,ds0,ds1,w) {
  n=length(e0)
  AICc0<-n*log(sum(e0^2)/n)+n*log(2*pi)+n*(n+sum(ds0))/(n-2-sum(ds0))
  wi<-as.matrix(w)
  e=t(e0*(1-wi)+wi*e1)
  ds=t(ds0*(1-wi)+wi*ds1)
  AICc0_1_i<-n*log(rowSums(e^2)/n)+n*log(2*pi)+n*(n+rowSums(ds))/(n-2-rowSums(ds))
  AICc0_1_i<1*AICc0

  # ff<-function(i){
  # wi<-as.matrix(w[i,])
  # #wi<-sqrt(wi)
  # #wi<-wi/sum(wi)
  # e=e0*(1-wi)+wi*e1
  # ds=ds0*(1-wi)+wi*ds1
  # AICc0_1_i<-n*log(sum(e^2)/n)+n*log(2*pi)+n*(n+sum(ds))/(n-2-sum(ds))
  # AICc0_1_i#<1*AICc0
  # }
  # sapply(1:n,ff)
}
