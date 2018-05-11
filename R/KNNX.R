#' KNNX
#' to be documented
#' @usage KNNX(X,Z=NULL,k_neighbors=NULL,h_bandwidth=NULL,diagnull=TRUE,
#' kernel='gauss_adapt',query=NULL,Zquery=NULL,ZC1=NULL,ZC2=NULL,
#' conditionZ=NULL,hd=0,norm=TRUE)
#' @param X  to be documented
#' @param Z  to be documented
#' @param k_neighbors  to be documented
#' @param h_bandwidth  to be documented
#' @param diagnull  to be documented
#' @param kernel  to be documented
#' @param query  to be documented
#' @param Zquery  to be documented
#' @param ZC1  to be documented
#' @param ZC2  to be documented
#' @param conditionZ  to be documented
#' @param hd  to be documented
#' @param norm  to be documented
#' @keywords internal
#' @return to be documented
KNNX<-function(X,Z=NULL,k_neighbors=NULL,h_bandwidth=NULL,diagnull=TRUE,kernel='gauss_adapt',query=NULL,Zquery=NULL,ZC1=NULL,ZC2=NULL,conditionZ=NULL,hd=0,norm=TRUE){
  X=as.matrix(X)
  if(is.null(query)) query=X
  if(!is.null(dim(X))) nr=nrow(X) else nr=length(X)
  k=k_neighbors+1;h=h_bandwidth;
  nb1<-knn(X,query,k = k)
  if(!is.null(dim(query))) m=nrow(query) else m=length(query)
  jj=as.numeric(t(nb1$nn.idx))
  ii=rep(seq_along(rep(k, m)), rep(k, m))
  long=length(ii)
  if(length(kernel)==1){
    if(!is.null(Z)){
      if(kernel=='bisq_Z') dd = bisq(Z[ii]-Zquery[jj],h) else if(kernel=='gauss_Z') dd = gauss_adapt(Z[ii]-Zquery[jj],h) else if(kernel=='bin_Z') dd = as.numeric(Z[ii]-Zquery[jj]<h)
    } else if(
      kernel=='knn') dd=rep(1/k,long)
    else if(kernel %in% c('bisq_knn','gauss_adapt','sheppard')){
      d <- as.numeric(t(nb1$nn.dists))
      DMmax=nb1$nn.dists[,k]*1.01
      dmax <-rep(DMmax,each=k)
      dmax[dmax==0]<-1
      if(kernel=='bisq_knn') dd = bisq(d,dmax) else if(kernel=='gauss_adapt') dd = gauss_adapt(d, dmax) else if(kernel=='sheppard') dd=(d - dmax)^2/(d * dmax) else dd=rep(1/k,long)
    }
    else if(kernel %in% c('gauss','bisq','bin')){
      if(kernel=='bisq') dd = bisq(d,h) else if(kernel=='gauss') dd = gauss_adapt(d,h) else if(kernel=='bin') dd = as.numeric(d<h)
    }
    if(!is.null(conditionZ)){
      if(is.null(ZC2)) ZC2=ZC1
      if(conditionZ=='=') {conditionZ=which(ZC1[ii]==ZC2[jj])}
      else if(conditionZ=='>') {conditionZ=which(ZC1[ii]>ZC2[jj])}
      else if(conditionZ=='<') {conditionZ=which(ZC1[ii]<ZC2[jj])}
      else if(conditionZ=='>=') {conditionZ=which(ZC1[ii]>=ZC2[jj])}
      else if(conditionZ=='<=') {conditionZ=which(ZC1[ii]<=ZC2[jj])}
      dd[conditionZ]<-dd[conditionZ]*hd
    }
    W <- sparseMatrix(i = ii, j =jj,  dims = c(m,nr), x = dd)
    if(diagnull) diag(W)=0
    if(norm) W<-normW(W)
  } else {
    if  (kernel[1]=='knn') dd=rep(1/k,long)
    else if(kernel[1] %in% c('bisq_knn','gauss_adapt')){
      d <- as.numeric(t(nb1$nn.dists))
      DMmax=nb1$nn.dists[,k]*1.01
      dmax <-rep(DMmax,each=k)
      dmax[dmax==0]<-1
      if(kernel[1]=='bisq_knn') dd = bisq(d,dmax) else if(kernel[1]=='gauss_adapt') dd = gauss_adapt(d, dmax) else if(kernel[1]=='sheppard') dd=(d - dmax)^2/(d * dmax) else dd=rep(1/k,long)
    }
    else if(kernel %in% c('gauss','bisq','bin')){
      if(kernel=='bisq') dd = bisq(d,h[1]) else if(kernel=='gauss') dd = gauss_adapt(d,h[1]) else if(kernel=='bin') dd = as.numeric(d<h[1])
    }
    if(!is.null(conditionZ)){
      if(is.null(ZC2)) ZC2=ZC1
      if(conditionZ=='=') {conditionZ=which(ZC1[ii]==ZC2[jj])}
      else if(conditionZ=='>') {conditionZ=which(ZC1[ii]>ZC2[jj])}
      else if(conditionZ=='<') {conditionZ=which(ZC1[ii]<ZC2[jj])}
      else if(conditionZ=='>=') {conditionZ=which(ZC1[ii]>=ZC2[jj])}
      else if(conditionZ=='<=') {conditionZ=which(ZC1[ii]<=ZC2[jj])}
      dd[conditionZ]<-dd[conditionZ]*hd
    }
    h2=h[length(h)]
    if(kernel[2]=='bisq_Z') dd2 = bisq(Z[ii]-Zquery[jj],h2) else if(kernel[2]=='gauss_Z') dd2 = gauss_adapt(Z[ii]-Zquery[jj],h2) else if(kernel[2]=='bin_Z') dd2 = as.numeric(Z[ii]-Zquery[jj]<h2)

    W <- sparseMatrix(i = ii, j =jj,  dims = c(m,nr), x = dd)
    W2 <- sparseMatrix(i = ii, j =jj,  dims = c(m,nr), x = dd2)
    W<-W*W2
    if(diagnull) diag(W)=0
    if(norm) W<-normW(W)
  }
  W
}
