#' mgwr_stage1
#' to be documented
#' @usage mgwr_stage1(Y,XV,XC,ALL_X,W,indexG=indexG,Wd,NN,isgcv,
#' TP,SE,Model,doMC,ncore)
#' @param Y  A vector
#' @param XV A matrix with covariates with spatially varying parameters
#' @param XC A matrix with covariates with stationnary parameters
#' @param indexG Precomputed Matrix of indexes of NN neighbors.
#' @param Wd Precomputed Matrix of weights.
#' @param NN Number of spatial Neighbours for kernels computations
#' @param isgcv leave one out cross validation, default FALSE
#' @param TP index of target points, default 1:n
#' @param SE If standard error are computed, default FALSE
#' @param Model character containing the type of model:
#'  Possible values are "OLS", "SAR", "GWR" (default), "MGWR" ,
#'   "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv", "MGWRSAR_0_kc_kv",
#'   "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' @param KernelTP  Kernel type for extrapolation of Beta from Beta(TP)
#' @param kWtp  Number of neighbours for extrapolation of Beta from Beta(TP)
#' @noRd
mgwr_stage1<-function(Y,XV,XC,ALL_X,W,indexG=indexG,Wd,NN,isgcv,TP,SE,Model,doMC,ncore){
  n=length(Y)
  ntp=length(TP)
  if(!is.null(XV)) m=ncol(XV) else m=0
  SY<- matrix(0,nrow=n, ncol=1)
  XX<-matrix(0,ncol=ncol(XC),nrow=n)
  if(isgcv) loo=-1 else loo=1:NN
  z=1

  if(doMC) { # & Sys.info()["sysname"]=='Windows'
    registerDoParallel(cores=ncore)
  } else registerDoSEQ()

  if(ncore>1) myblocks<-split(1:length(TP), ceiling(seq_along(TP)/round(length(TP)/ncore))) else myblocks<-list(b1=1:length(TP))

  #res<-foreach(z =1:length(TP),.combine="comb",.inorder=FALSE)  %dopar% {
  res<-foreach(myblock =1:length(myblocks),.combine="comb",.inorder=FALSE)  %dopar% {
    for(z in myblocks[[myblock]]){
      i=TP[z]
      index=indexG[z,loo]
      wd<-sqrt(Wd[z,loo])
      Yw<-wd*Y[index]
      if(Model!='MGWRSAR_1_kc_0') Xw=wd*XV[index,] else {
        PhWy<-try(PhWY_C(as.matrix(Y[index]),as.matrix(ALL_X[index,]*wd),W[index,index],rep(1,length(index))),silent =TRUE)
        if(is(PhWy,'try-error')) PhWy=PhWY_R(as.matrix(Y[index]),as.matrix(ALL_X[index,]*wd),W[index,index],rep(1,length(index)))
        Xw=PhWy*wd
        XV<-W%*%Y
        XVi<-XV[i]
      }
      if(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv')){
        PhWy=PhWY_R(as.matrix(Y[index]),as.matrix(ALL_X[index,]*wd),W[index,index],rep(1,length(index)))
        Xw=cbind(Xw,PhWy*wd)
        XVi<-c(XV[i,],(W%*%Y)[i])
      } else XVi<-XV[i,]
      matB<-QRcpp2_C(Xw,as.matrix(wd*Y[index]),as.matrix(wd*XC[index,]))
      xx<- XVi %*% matB$XCw
      sy<-as.numeric(XVi %*% matB$SY)
      if(z==myblocks[[myblock]][1]){
        xxz=matrix(xx,nrow=1)
        syz=matrix(sy,nrow=1)
      } else {
        xxz=rbind(xxz,matrix(xx,nrow=1))
        syz=rbind(syz,matrix(sy,nrow=1))
      }
    }
    rm(index,wd,Yw,Xw,matB,XVi,xx,sy)
    gc()
    list(xx=xxz,sy=syz)
  }
  XX[TP,]=res$xx
SY[TP]=res$sy
  if(ntp<length(Y)){
    # Wtp<- normW(Matrix::t(sparseMatrix(i = rep(1:ntp,each=NN), j = as.numeric(t(indexG)),  dims = c(ntp,n), x =as.numeric(t(Wd))))[-TP,]) ## revoir
    #XX[-TP,]=as.matrix(Wtp%*% XX[TP,])
    #SY[-TP]=as.numeric(Wtp%*% SY[TP])
    SY=Y[TP]-SY[TP]
    XX=XC[TP,]-XX[TP,]
  } else {
    SY=Y-SY
    XX=XC-XX
  }
  model1<-lm(SY~XX-1)
  Betac<-coefficients(model1)
  if(any(is.na(Betac))){
    cat(paste0('Warnings: coefficients of ',names(Betac)[is.na(Betac)],' has been set to 0 du to collinearity'))
    Betac[is.na(Betac)]<-0
  }
  if(SE) se <- sqrt(diag(vcov(model1))) else se=NULL
  ZZ<-Y-XC %*%as.matrix(Betac)
  list(ZZ=ZZ,Betac=Betac,se=se)
}
