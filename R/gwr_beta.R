#' gwr_beta
#' to be documented
#' @usage gwr_beta(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,k
#' ernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,
#' TP_estim_as_extrapol=FALSE,get_ts=FALSE,get_s=FALSE,
#' get_Rk=FALSE,TP_cor=NULL)
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
#' @param ncore  Number of cores for parallel computation.
#' @param TP_estim_as_extrapol  Boolean for prediction mode.
#' @param get_ts Boolean for computing Trace(S)
#' @param get_Rk Boolean for computing Trace(S)
#' @noRd
#' @return A list with Betav, standard error, edf and trace(hatMatrix)
gwr_beta<-function(Y,XV,ALL_X,TP,indexG,Wd,NN,W=NULL,isgcv=FALSE,SE=FALSE,kernels=NULL,H=NULL,adaptive=NULL,doMC=FALSE,ncore=1,TP_estim_as_extrapol=FALSE,get_ts=FALSE,get_s=FALSE,get_Rk=FALSE,TP_cor=NULL)
{
  if(get_s) get_ts=TRUE
  Rk<-Rkk<-list()
  n=length(Y)
  if(!is.null(TP_cor)) {
    TP<-TP_cor
    some_null_beta<-TRUE
  } else some_null_beta<-FALSE
  ntp=length(TP)
  if(TP_estim_as_extrapol) {
    SE=FALSE
    isgcv=FALSE
  }
  if(get_ts | get_s | SE) {
    tS<-0
    TS=rep(0,n)
  }
  if(get_Rk){ Rk<-array(0,dim = c(n,n,ncol(XV)),dimnames=list(NULL,NULL,colnames(XV)))
  }
  if(get_s)  Shat <- matrix(0,ncol=n,nrow=n) else Shat=NULL # a corriger
  if(get_s | SE)  SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), ncol(XV), ncol(XV) + 1)) else SEV=NULL # a corriger
  if(!is.null(XV)) m=ncol(XV) else m=0
  tS<-0
  namesXV=colnames(XV)
  if (!is.null(W)) {
    PhWy=PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1,n))
    XV = cbind(XV,PhWy)
  }
  if(isgcv) loo=-1 else loo=1:NN
  if(doMC) {
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()
  if(ncore>1) myblocks<-split(1:length(TP), ceiling(seq_along(TP)/round(length(TP)/ncore))) else myblocks<-list(b1=1:length(TP)) ## myblocks and z index over 1:length(TP)
  edf=n;

  if(TP_estim_as_extrapol) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))


  if(length(myblocks[[length(myblocks)]])==1){
    myblocks[[length(myblocks)-1]]<-c(myblocks[[length(myblocks)-1]],myblocks[[length(myblocks)]])
    myblocks=myblocks[-length(myblocks)]
  }
  res<-foreach(myblock =1:length(myblocks),.combine="comb",.inorder=FALSE)  %dopar% {
    ts=c()
    for(z in myblocks[[myblock]]){
      #if(z==3) browser()
      #cat("z =",z,' ')
      correc_lambda=FALSE
      index=indexG[z,loo]
      wd<-sqrt(Wd[z,loo])
      Yw<-wd*Y[index]
      Xw=as.matrix(wd*XV[index,])
      lml=lm.fit(as.matrix(Xw),as.matrix(Yw))
      betav=lml$coefficients
      coefNA<-which(is.na(betav))
      betav[coefNA]<-0
      if(length(coefNA)>0 & length(coefNA)<ncol(Xw)) {
        lml=lm.fit(as.matrix(Xw[,-coefNA]),as.matrix(Yw))
        betav[-coefNA]<-lml$coefficients
      }
      if(!is.null(W) & abs(betav[m + 1])>1) {
        correc_lambda=TRUE
        betav[m + 1]=sign(betav[m + 1])*0.99
        if(m>0) {
          lml=lm.fit(as.matrix(Xw[,-c(coefNA,m+1)]),as.matrix(Yw-betav[m + 1]*Xw[,m + 1]))
          betav[setdiff(1:m,coefNA)]=lml$coefficients
        }
      }
      if(get_ts) {
        tS=tS+lm.influence(lml)$hat[1]
        ts=c(ts,lm.influence(lml)$hat[1])
      }
      if(get_s) {
        coef_NON_NA=setdiff(1:ncol(Xw),coefNA)
        if(length(coef_NON_NA)>0){
        XwX<-try(solve(crossprod(Xw[, coef_NON_NA], Xw[, coef_NON_NA]),silent = TRUE))
        Zwi =try( XwX %*% t(Xw[, coef_NON_NA]),silent = TRUE)
        Shat[TP[z],index]<-((XV[TP[z], coef_NON_NA] %*% Zwi)*wd)
        if(get_Rk) {
          for(nx in coef_NON_NA) Rk[TP[z],index,nx]<- ((XV[TP[z], nx] * Zwi[nx, ])*wd)
        }
        } else  Shat[TP[z],index]<-0
      }
      if(SE & !isgcv) {
        coef_NON_NA=setdiff(1:ncol(Xw),coefNA)
        if(correc_lambda) SEV[TP[z],coef_NON_NA]<- c(sqrt(dvcov(lml)),NA) else SEV[TP[z],coef_NON_NA]<- sqrt(dvcov(lml))
        ##if(correc_lambda) SEV[z,coef_NON_NA] <- c(sqrt(dvcov(lml)),NA) else SEV[z,coef_NON_NA]<- sqrt(dvcov(lml))
      }
      #if(!TP_estim_as_extrapol){ Betav[TP[z],]<-betav } else {Betav[z,]<-betav}
      #Betav[z,]<-betav

      Betav[TP[z],]<-betav # si TP= 1:n, TP[z]=z ; si TP
    }
    # if(!(SE & !isgcv)) {
    #   sev=NULL
    # } else {
    #   sev=SEV[TP[myblocks[[myblock]]],] ## a corriger
    # }
    if(get_s)  Shat=Shat[TP[myblocks[[myblock]]],] else Shat=NULL
    if(get_Rk)  {
      for(nx in 1:ncol(XV)) {
        Rkk[[nx]]<-Rk[TP[myblocks[[myblock]]],,nx]
      }
      }

    if(length(ts)>0) TS=as.matrix(ts,ncol=1) else TS=NULL
    rm(index,wd,Yw,Xw,betav)
    if(SE) sev<-as.matrix(SEV[TP[myblocks[[myblock]]],],ncol=m) else sev=NULL
    list(betav=as.matrix(Betav[TP[myblocks[[myblock]]],],ncol=m),sev=sev,tS=tS,Shat=Shat,TS=TS,Rk=Rkk) ## return foreach
  }
  if(m==1) res$betav=as.numeric(t(res$betav))
  if(TP_estim_as_extrapol) Betav=matrix(0,nrow=ntp,ncol= ifelse(is.null(W), m, m + 1)) else Betav=matrix(0,nrow=n,ncol= ifelse(is.null(W), m, m + 1))
  if(!TP_estim_as_extrapol) {
    Betav[TP,]<-res$betav
    if(SE) {
      edf=n;
      SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W), m, m + 1))
      SEV[TP,]=res$sev
    }
    if(get_ts)  {
      tS=sum(res$tS)
      TS[TP]=as.numeric(res$TS)
      #names(TS)=TP
    }
    if(get_s)  {
      Shat=res$Shat
      if(get_Rk) {
        Rk=list()
        for(nx in 1:ncol(XV)){
          index=1:length(myblocks)+(nx-1)*length(myblocks)
          rk=res$Rk[[index[1]]]
          for(zz in index[-1]){
            rk<-rbind(rk,res$Rk[[zz]])
          }
          Rk[[colnames(XV)[nx]]]<-rk
        }
      }
    }
  } else {
    Betav<-res$betav
  }
  if(SE) colnames(SEV)=colnames(XV)
  # if(ntp<length(Y) & !TP_estim_as_extrapol & !some_null_beta){
  #   #browser()
  #   Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=n), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,])
  #   Betav[-TP,]=as.matrix(Wtp%*% Betav[TP,])
  #   if(SE) SEV[-TP,]=as.matrix(Wtp%*% SEV[TP,])
  # }
  if(is.null(W))  colnames(Betav)=namesXV else colnames(Betav)=c(namesXV,'lambda')

  if(get_s | get_ts | SE) list(Betav=Betav,SEV=SEV,edf=n-tS,tS=tS,Shat=Shat,TS=TS,Rk=Rk) else list(Betav=Betav,SEV=NULL,edf=NULL,tS=NULL,Shat=NULL,TS=NULL)
}

