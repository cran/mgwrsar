#' PhWY_R
#' first stage of S2SLS in case
#' @usage PhWY_R(Y, X,W, Wi=rep(1,nrow(X)))
#' @param Y  A vector of responses
#' @param X  A matrix of covariates
#' @param W  A weight matrix for spatial autocorrelation
#' @param Wi A vectof weights
#' @noRd
#' @return to be documented
PhWY_R<-function(Y, X,W, Wi=rep(1,nrow(X))){
Wy=as.numeric(W%*%Y)
WX<-W%*%X[,-1]
WWX<-W%*%WX
WWWX<-W%*%WWX
Z<-as.matrix(cbind(X,WX,WWX,WWWX))
return(fitted(lm.fit(Z,Wy)))
}

