#' pastlm
#' to be documented
#' @usage pastlm(formula,data,K,doMC=TRUE,ncore=8)
#' @param formula to be documented
#' @param data to be documented
#' @param K to be documented
#' @param doMC to be documented
#' @param ncore to be documented
#' @noRd
#' @return a vector of weights.
pastlm<-function(formula,data,K,doMC=TRUE,ncore=8){
  ## data must be arrange by time
  idx=(K+2):nrow(data)
  Beta=matrix(NA,nrow=nrow(data),ncol=K)
  if(doMC) {
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()
  beta<-foreach(x =1:length(idx),.combine="rbind")  %dopar% {
     id=1:idx[x]
     coef(lm(formula,data[id,]))
   }
  beta[is.na(beta)]<-0
  Beta[idx,]<-beta
  Beta
}
