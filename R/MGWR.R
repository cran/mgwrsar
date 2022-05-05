#' MGWR
#' to be documented
#' @usage MGWR(Y,XC,XV,ALL_X=NULL,S,H,NN, kernels,adaptive=F,Type="GD",
#' SE=FALSE,isgcv=F,W=NULL,remove_local_outlier=FALSE,outv=0,TP=NULL,
#' KernelTP='sheppard',kWtp=8,Model,indexG=NULL,Wd=NULL,dists=NULL,
#' doMC=FALSE,ncore=1,S_out=NULL)
#' @param Y  A vector
#' @param XC  A matrix with covariates with stationnary parameters
#' @param XV   A matrix with covariates with spatially varying parameters
#' @param S  A matrix with variables used in kernel
#' @param H  A vector of bandwidths
#' @param NN Number of spatial Neighbours for kernels computations.
#' @param kernels  A vector of kernel types
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel.
#' @param Type  Type of Genelarized kernel product ('GD' only
#'  spatial,'GDC' spatial + a categorical variable,
#'  'GDX' spatial + a continuous variable and other
#'   combinations like 'GDXXC','GDXCC',...)
#' @param SE  If standard error are computed, default FALSE
#' @param isgcv  leave one out cross validation, default FALSE
#' @param W  A weight matrix for spatial autocorrelation
#' @param remove_local_outlier Remove local outlier
#' @param outv  A treshold for removing local outlier
#' @param TP  index of target points, default NULL
#' @param indexG  Precomputed Matrix of indexes of NN neighbors, default NULL.
#' @param dists  Precomputed Matrix of spatial distances, default NULL
#' @param Model character containing the type of model:
#'  Possible values are "OLS", "SAR", "GWR" (default), "MGWR" ,
#'   "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv", "MGWRSAR_0_kc_kv",
#'   "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' @return a list of object for MGWRSAR wrapper
#' @noRd
MGWR<-function(Y,XC,XV,ALL_X=NULL,S,H,NN, kernels,adaptive=F, Type = "GD",SE=FALSE, isgcv=F,W=NULL,remove_local_outlier=FALSE,outv=0,TP=NULL,Model,indexG=NULL,Wd=NULL,dists=NULL,doMC=FALSE,ncore=1,S_out=NULL){
  se = NULL
  sev = NULL
  coord=S[,1:2]
  if(ncol(S)>2) Z=matrix(S[,3:ncol(S)]) else Z=NULL
  if (!is.null(XC)) XC <- as.matrix(XC)
  if (!is.null(XV)) XV <- as.matrix(XV)
  ALL_X = cbind(XC, XV)
  n <- NROW(Y)
  m <- ncol(XV)
  K <- ncol(XC)
  if (Model %in% c("MGWRSAR_0_kc_kv", "MGWRSAR_0_0_kv")) {
   PhWy=PhWY_R(as.matrix(Y), as.matrix(ALL_X), W, rep(1,n))
    if (Model == "MGWRSAR_0_kc_kv")
      XC = cbind(XC, PhWy)
    else XC = as.matrix(PhWy,ncol=1)
  }
  ### prep wd
  ########
  if(is.null(Wd)){
    if(is.null(S_out)) {
      Z=S[TP,]
      pred=FALSE
    } else {
      Z=S_out
      pred=TRUE
    }
      stage1=prep_w(H=H,kernels=kernels,coord_i=Z,coord_j=S,NN=NN,ncolX=ncol(XV),Type=Type,adaptive=adaptive,dists=dists,indexG=indexG,rowNorm=TRUE)
      indexG=stage1$indexG
      dists=stage1$dists
      Wd=stage1$Wd
  }

  mgwr1=mgwr_stage1(Y=Y,XV=XV,XC=XC,ALL_X=ALL_X,W=W,indexG=indexG,Wd=Wd,NN=NN,isgcv=isgcv,TP=TP,SE=SE,Model=Model,doMC=doMC,ncore=ncore)

  if (Model %in% c("MGWRSAR_1_kc_kv","MGWRSAR_1_kc_0")) {
    W2=W
  } else W2=NULL
  model=gwr_beta(Y=mgwr1$ZZ,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,W=W2,isgcv=isgcv,SE=SE,remove_local_outlier=remove_local_outlier,outv=outv,doMC=doMC,ncore=ncore, pred=pred)

  if(SE & !isgcv) list(Betac=mgwr1$Betac,Betav=model$Betav,SEV=model$SEV,se=mgwr1$se,edf=n-model$tS-length(mgwr1$Betac),tS=model$tS+length(mgwr1$Betac)) else list(Betac=mgwr1$Betac,Betav=model$Betav,SEV=NULL,edf=NULL,tS=NULL)
}
