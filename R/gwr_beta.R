#' gwr_beta
#' to be documented
#' @usage gwr_beta(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,
#' KernelTP='sheppard',kWtp=8,doMC=FALSE,ncore=1,pred=FALSE,get_ts=FALSE)
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
#' @noRd
#' @return A list with Betav, standard error, edf and trace(hatMatrix)
gwr_beta<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,pred=FALSE,get_ts=FALSE)
{
  n=length(Y)
  ntp=length(TP)
  if(pred) {
    SE=FALSE
    isgcv=FALSE
  }
  if(!is.null(XV)) m=ncol(XV) else m=0
  tS=0
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
    if(SE) {
      tS<-0;
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
    }
    if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
    for(z in myblocks[[myblock]]){
      index=indexG[z,loo]
      wd<-sqrt(Wd[z,loo])
      Yw<-wd*Y[index]
      Xw=as.matrix(wd*XV[index,])
      lml=lm.fit(as.matrix(Xw),as.matrix(Yw))
      betav=lml$coefficients
      coefNA<-which(is.na(betav))
      betav[coefNA]<-0
      if(length(coefNA)>0) {
        lml=lm.fit(as.matrix(Xw[,-coefNA]),as.matrix(Yw))
        betav[-coefNA]<-lml$coefficients
      }
      if(!is.null(W) & abs(betav[m + 1])>1) {
        betav[m + 1]=sign(betav[m + 1])*0.99
        if(m>0) {
          lml=lm.fit(as.matrix(Xw[,-c(coefNA,m+1)]),as.matrix(Yw-betav[m + 1]*Xw[,m + 1]))
          betav[setdiff(1:m,coefNA)]=lml$coefficients
        }
      }
      if(get_ts & !SE){
        coef_NON_NA=setdiff(1:ncol(Xw),coefNA)
        Zwi = try(solve(crossprod(Xw[, coef_NON_NA], Xw[, coef_NON_NA])) %*% t(Xw[, coef_NON_NA]),silent = TRUE)
        tS=tS+ifelse(class(Zwi)[1]=='try-error',0,(Xw[1, coef_NON_NA] %*% Zwi)[, 1])
      }
      if(SE & !isgcv) {
        coef_NON_NA=setdiff(1:ncol(Xw),coefNA)
        rss <- sum(lml$residuals^2)
        rdf <- length(Yw) - ncol(Xw)+length(coefNA)
        resvar <- rss/rdf
        R <- chol2inv(lml$qr$qr)
        diagR=diag(R)
        SEV[TP[z],coef_NON_NA] <- sqrt(diagR * resvar)
        Zwi = try(solve(crossprod(Xw[, coef_NON_NA], Xw[, coef_NON_NA])) %*% t(Xw[, coef_NON_NA]),silent = TRUE)
        tS=tS+ifelse(class(Zwi)[1]=='try-error',0,(Xw[1, coef_NON_NA] %*% Zwi)[, 1])
      }
      if(!pred){ Betav[TP[z],]<-betav } else {Betav[z,]<-betav}
    }
    if(!(SE & !isgcv)) {
      sev=NULL
      if(!get_ts) tS=NULL
    } else {
      sev=SEV[TP[myblocks[[myblock]]],]
    }
    rm(index,wd,Yw,Xw,betav)
    gc()
    list(betav=Betav[TP[myblocks[[myblock]]],],sev=sev,tS=tS)
  }

  if(pred) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))

  if(!pred) {
    Betav[TP,]<-res$betav
    if(SE) {
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
      SEV[TP,]=res$sev
      tS=sum(res$tS)
    }
    if(get_ts)  tS=sum(res$tS)
  } else {
    Betav<-res$betav
  }
  if(SE) colnames(SEV)=colnames(XV)
  if(ntp<length(Y) & !pred){
    Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=NN), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,])
    Betav[-TP,]=as.matrix(Wtp%*% Betav[TP,])
    if(SE) SEV[-TP,]=as.matrix(Wtp%*% SEV[TP,])
  }
  if(is.null(W))  colnames(Betav)=namesXV else colnames(Betav)=c(namesXV,'lambda')
  if(SE & !isgcv & !pred) list(Betav=Betav,SEV=SEV,edf=n-tS,tS=tS)  else if(get_ts) list(Betav=Betav,SEV=NULL,edf=NULL,tS=tS) else list(Betav=Betav,SEV=NULL,edf=NULL,tS=NULL)
}

