#' gwr_beta
#' to be documented
#' @usage gwr_beta(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=F,SE=FALSE,
#' remove_local_outlier=F,outv=0.01,KernelTP='sheppard',kWtp=8,
#' doMC=FALSE,ncore=1)
#' @param Y A vector of response
#' @param XV A matrix with covariates with non stationnary parameters
#' @param ALL_X A matrix with all covariates
#' @param coord A matrix with spatial coordinates
#' @param TP An index of target points.
#' @param indexG Precomputed Matrix of indexes of NN neighbors.
#' @param Wd Precomputed Matrix of weights.
#' @param NN Number of spatial Neighbours for kernels computations
#' @param W The spatial weight matrix for spatial dependence
#' @param isgcv leave one out cross validation, default FALSE
#' @param SE If standard error are computed, default FALSE
#' @param remove_local_outlier  Remove local outlier
#' @param outv percentile threshold for outlier.
#' @param KernelTP  Kernel type for extrapolation of Beta from Beta(TP)
#' @param kWtp  Number of neighbours for extrapolation of Beta from Beta(TP)
#' @param doMC  Boolean for parallel computation.
#' @param pred  Is the GWR used for prediction for target points ?
#' @param mstop   Number of iterations for mboost.
#' @param nu  Learning rate for mboost.
#' @noRd
#' @return A list with Betav, standard error, edf and trace(hatMatrix)
gwr_beta_glmboost1<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=F,SE=FALSE,remove_local_outlier=F,outv=0.01,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,pred=FALSE,mstop=150,nu=0.1)
{
  n=length(Y)
  ntp=length(TP)
  if(pred) {
    SE=FALSE
    isgcv=F
  }
  if(!is.null(XV)) m=ncol(XV) else m=0
  tS=0
  namesXV=colnames(XV)
  if (!is.null(W)) {
    PhWy=PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1,n))
    XV = cbind(XV,PhWy)
  }## ce point pose probleme dans le cas isgcv

  if(isgcv) loo=-1 else loo=1:NN
  if(doMC) {
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()
  if(ncore>1) myblocks<-split(1:length(TP), ceiling(seq_along(TP)/round(length(TP)/ncore))) else myblocks<-list(b1=1:length(TP))
  res<-foreach(myblock =1:length(myblocks),.combine="comb",.inorder=FALSE)  %dopar% {
    if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
    for(z in myblocks[[myblock]]){
      index=indexG[z,loo] #### commencer ici les adaptations boost
      if(remove_local_outlier | !is.null(W) | SE){
        stop('remove_local_outlier not implement with boosting adn SAR and SE')
      }
      betav<-rep(0,ncol(XV))
      names(betav)<-colnames(XV)
      ## version1
      res=glmboost(x=as.matrix(XV[index,]), y=as.numeric(Y[index]),weights=Wd[z,loo],center=TRUE,control = boost_control(mstop = mstop,nu=nu))
      mycoef<-coef(res,off2int = T)
      betav[names(mycoef)]<-mycoef
      if(!pred){ Betav[TP[z],]<-betav} else {Betav[z,]<-betav}
    }
    if(!(SE & !isgcv)) {
      sev=NULL
      tS=NULL
    } else {
      sev=SEV[TP[myblocks[[myblock]]],]
    }
    rm(index,betav)
    gc()
    list(betav=Betav[TP[myblocks[[myblock]]],],sev=sev,tS=tS)
  } #

  if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))

  if(!pred) {
    Betav[TP,]<-res$betav
  } else {
    Betav<-res$betav
  }

  if(SE) colnames(SEV)=colnames(XV)
  if(ntp<length(Y) & !pred){
    #if(KernelTP=='Wd'){ ### on ne recalcule pas W
    Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=NN), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,])
    Betav[-TP,]=as.matrix(Wtp%*% Betav[TP,])
    if(SE) SEV[-TP,]=as.matrix(Wtp%*% SEV[TP,])
  }
  if(is.null(W))  colnames(Betav)=namesXV else colnames(Betav)=c(namesXV,'lambda')
  if(SE & !isgcv & !pred) list(Betav=Betav,SEV=SEV,edf=n-tS,tS=tS) else list(Betav=Betav,SEV=NULL,edf=NULL,tS=NULL)
}



