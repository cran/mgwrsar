#' gwr_beta_glm
#' to be documented
#' @usage gwr_beta_glm(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,
#' KernelTP='sheppard',kWtp=8,doMC=FALSE,ncore=1,pred=FALSE,get_ts=FALSE
#' family=NULL)
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
#' @param kWtp  Number of neighbours for extrapolation of Beta from Beta(TP)
#' @param doMC  Boolean for parallel computation.
#' @param ncore  Number of cores for parallel computation.
#' @param pred  Boolean for storing fitted value.
#' @param get_ts Boolean for computing Trace(S)
#' @param family 	a Family object see(glmboost help)
#' @noRd
#' @return A list with Betav, standard error, edf and trace(hatMatrix)
gwr_beta_glm<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,pred=FALSE,get_ts=FALSE,family=NULL)
{
  if(is.null(family)) family = gaussian(link = "identity")
  n=length(Y)
  ntp=length(TP)
  if(pred) {
    SE=FALSE
    isgcv=FALSE
  }
  if(!is.null(XV)) m=ncol(XV) else m=0
  tS=0
  namesXV=colnames(XV)
  if(isgcv) loo=-1 else loo=1:NN
  if(doMC) {
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()
  if(ncore>1) myblocks<-split(1:length(TP), ceiling(seq_along(TP)/round(length(TP)/ncore))) else myblocks<-list(b1=1:length(TP))
  res<-foreach(myblock =1:length(myblocks),.combine="comb",.inorder=FALSE)  %dopar% {
    if(SE) {
      tS<-0;
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
    }
    if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
    for(z in myblocks[[myblock]]){
      index=indexG[z,loo]
      wd<-Wd[z,loo]
      dataglm<-data.frame(Y=as.matrix(Y[index]),as.matrix(XV[index,]))
      lml=glm(formula=as.formula('Y~.-1'),data=dataglm,family=family,weights=wd)
      betav=lml$coefficients
      coefNA<-which(is.na(betav))
      betav[coefNA]<-0
      if(SE & !isgcv) {
        coef_NON_NA=setdiff(1:ncol(XV),coefNA)
        SEV[TP[z],coef_NON_NA] <- sqrt(diag(vcov(lml)))
      }
      if(!pred){ Betav[TP[z],]<-betav } else {Betav[z,]<-betav}
    }
    if(!(SE & !isgcv)) {
      sev=NULL
    } else {
      sev=SEV[TP[myblocks[[myblock]]],]
    }
    rm(index,wd,betav)
    gc()
    list(betav=Betav[TP[myblocks[[myblock]]],],sev=sev,tS=NULL)
  }

  if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))

  if(!pred) {
    Betav[TP,]<-res$betav
    if(SE) {
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
      SEV[TP,]=res$sev
    }
  } else {
    Betav<-res$betav
  }
  if(SE) colnames(SEV)=colnames(XV)
  if(ntp<length(Y) & !pred){
    Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=NN), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,])
    Betav[-TP,]=as.matrix(Wtp%*% Betav[TP,])
    if(SE) SEV[-TP,]=as.matrix(Wtp%*% SEV[TP,])
  }
  colnames(Betav)=namesXV
  if(SE & !isgcv & !pred) list(Betav=Betav,SEV=SEV,edf=NULL,tS=NULL) else list(Betav=Betav,SEV=NULL,edf=NULL,tS=NULL)
}

