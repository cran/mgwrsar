#' Internal GWR Estimation via Boosting (mboost)
#'
#' @description
#' This function estimates GWR coefficients using Component-wise Gradient Boosting (`mboost`).
#' It is particularly useful for high-dimensional data or when variable selection is required
#' locally. It fits a `glmboost` model at each target location.
#'
#' @usage gwr_beta_glmboost(Y, XV, ALL_X, TP, indexG, Wd, NN, W = NULL, isgcv = FALSE,
#' SE = FALSE, kernels = NULL, H = NULL, adaptive = NULL, ncore = 1,
#' TP_estim_as_extrapol = FALSE, mstop = 150, nu = 0.1, family = NULL)
#'
#' @param Y A numeric vector of the response variable.
#' @param XV A matrix of covariates with spatially varying parameters.
#' @param ALL_X A matrix containing all covariates (used for SAR computations).
#' @param TP A vector of indices for target points.
#' @param indexG A precomputed matrix of indices of Nearest Neighbours.
#' @param Wd A precomputed matrix of spatial weights.
#' @param NN An integer indicating the number of spatial neighbours.
#' @param W A spatial weight matrix for spatial dependence (optional).
#' @param isgcv Logical; Not used for boosting (CV is usually internal to mboost), default FALSE.
#' @param SE Logical; Standard errors are NOT available for boosting methods. Must be FALSE.
#' @param kernels Character vector of kernel types (unused inside, kept for compatibility).
#' @param H Vector of bandwidths (unused inside, kept for compatibility).
#' @param adaptive Logical/Vector (unused inside, kept for compatibility).
#' @param ncore Integer; number of cores for parallel computation. Default is 1.
#' @param TP_estim_as_extrapol Logical; if TRUE, prediction mode is enabled.
#' @param mstop Integer; Number of boosting iterations (mboost parameter). Default is 150.
#' @param nu Numeric; Learning rate or step size (mboost parameter). Default is 0.1.
#' @param family An object of class `Family` from the `mboost` package (e.g., `Gaussian()`, `Binomial()`).
#'        If NULL, defaults to `Gaussian()`.
#'
#' @return A list containing:
#'   \item{Betav}{Matrix of local coefficients.}
#'   \item{SEV}{NULL (SE not supported for boosting).}
#'   \item{edf}{NULL.}
#'   \item{tS}{Fixed to 1 (Trace S not computed for boosting).}
#'
#' @seealso \code{\link[mboost]{glmboost}}, \code{\link[mboost]{boost_control}}
#' @keywords internal
#' @noRd
gwr_beta_glmboost<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,kernels=NULL,H=NULL,adaptive=NULL,ncore=1,TP_estim_as_extrapol=FALSE,mstop=150,nu=0.1,family=NULL)
{
  if(SE) stop('No variance estimation for glmboost and gamboost model')
  if(is.null(family) | family$family=="gaussian") family=Gaussian()
  n=length(Y)
  ntp=length(TP)
  if(!is.null(XV)) m=ncol(XV) else m=0
  namesXV=colnames(XV)
  if (!is.null(W)) {
    PhWy=PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1,n))
    XV = cbind(XV,PhWy)
  }

  if(isgcv) loo=-1 else loo=1:NN
  if(ncore>1) {
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()
  if(ncore>1) myblocks<-split(1:length(TP), ceiling(seq_along(TP)/round(length(TP)/ncore))) else myblocks<-list(b1=1:length(TP))
  res<-foreach(myblock =1:length(myblocks),.combine="comb",.inorder=FALSE)  %dopar% {
    if(TP_estim_as_extrapol) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
    for(z in myblocks[[myblock]]){
      index=indexG[z,loo]
      betav<-rep(0,ncol(XV))
      names(betav)<-colnames(XV)
      res=glmboost(x=as.matrix(XV[index,]), y=as.numeric(Y[index]),weights=Wd[z,loo],center=TRUE,control = boost_control(mstop = mstop,nu=nu),family = family)
      #cv5f <- cv(model.weights(res), type = "kfold",B=5)
      #cvm <- cvrisk(res, folds = cv5f, mc.cores = 1)
      #mstop0=max(10,mstop(cvm))
      #res=glmboost(x=as.matrix(XV[index,]), y=as.numeric(Y[index]),weights=Wd[z,loo],center=TRUE,control = boost_control(mstop = mstop0,nu=nu),family = family)
      #cat('i=',z,' mstop=',mstop0,' ')
      mycoef<-coef(res,off2int = T)
      names(mycoef)[names(mycoef)=='(Intercept)']<-'Intercept'
      betav[names(mycoef)]<-mycoef
      if(!TP_estim_as_extrapol){ Betav[TP[z],]<-betav} else {Betav[z,]<-betav}
    }
    #cat('\n\n')
    rm(index,betav)
    gc()
    list(betav=Betav[TP[myblocks[[myblock]]],])
  }
  if(TP_estim_as_extrapol) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
  if(!TP_estim_as_extrapol) {
    Betav[TP,]<-res$betav
  } else {
    Betav<-res$betav
  }
  if(ntp<length(Y) & !TP_estim_as_extrapol){
    Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=NN), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,])
    Betav[-TP,]=as.matrix(Wtp%*% Betav[TP,])
  }
  if(is.null(W))  colnames(Betav)=namesXV else colnames(Betav)=c(namesXV,'lambda')
  list(Betav=Betav,SEV=NULL,edf=NULL,tS=1)
}



