#' gwr_beta_glmboost
#' to be documented
#' @usage gwr_beta_glmboost((Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,
#' SE=FALSE,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,
#' TP_estim_as_extrapol=FALSE, mstop=150,nu=0.1,family=NULL)
#' @param Y A vector of response
#' @param XV A matrix with covariates with non stationnary parameters
#' @param ALL_X A matrix with all covariates
#' @param TP An index of target points.
#' @param indexG Precomputed Matrix of indexes of NN neighbors.
#' @param Wd Precomputed Matrix of weights.
#' @param NN Number of spatial Neighbours for kernels computations
#' @param W The spatial weight matrix for spatial dependence
#' @param isgcv leave one out cross validation, default FALSE
#' @param SE If standard error are computed, default FALSE
#' @param KernelTP  Kernel type for extrapolation of Beta from Beta(TP)
#' @param doMC  Boolean for parallel computation.
#' @param TP_estim_as_extrapol  Is the GWR used for prediction for
#' target points ?
#' @param mstop   Number of iterations for mboost.
#' @param nu  Learning rate for mboost.
#' @param family 	a Family object see(glmboost help)
#' @noRd
#' @return A list with Betav, standard error, edf and trace(hatMatrix)
gwr_beta_glmboost<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,TP_estim_as_extrapol=FALSE,mstop=150,nu=0.1,family=NULL)
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
  if(doMC) {
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



