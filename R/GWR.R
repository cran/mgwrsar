#' GWR
#' to be documented
#' @usage GWR(Y,XV,ALL_X,S,H,NN, kernels,adaptive=FALSE, Type = "GD",
#' SE=FALSE, isgcv=FALSE,W=NULL,TP=NULL,dists=NULL,indexG=NULL,
#' Wd=NULL,doMC=FALSE,ncore=1,Model,S_out=NULL,get_ts=FALSE)
#' @param Y  A vector of response
#' @param XV A matrix with covariates with stationnary parameters
#' @param ALL_X  A matrix with all covariates (XC,XV)
#' @param S  A matrix with variables used in kernel
#' @param H  A vector of bandwidths
#' @param NN Number of spatial Neighbours for kernels computations.
#' @param kernels  A vector of kernel types
#' @param adaptive  A vector of boolean to choose adaptive version for each kernel.
#' @param Type  Type of Genelarized kernel product ('GD' only spatial,'GDC' spatial
#' + a categorical variable,'GDX' spatial + a continuous variable,
#' and other combinations 'GDXXC','GDXCC',...)
#' @param SE  If standard error are computed, default FALSE
#' @param isgcv  leave one out cross validation, default FALSE
#' @param W  A weight matrix for spatial autocorrelation
#' @param TP  index of target points, default NULL
#' @param dists  Precomputed Matrix of spatial distances, default NULL
#' @param indexG  Precomputed Matrix of indexes of NN neighbors, default NULL.
#' @param Wd Precomputed Matrix of weights.
#' @param doMC  Boolean for parallel computation.
#' @param ncore  Number of cores for parallel computation.
#' @param Model character containing the type of model:
#'  Possible values are "OLS", "SAR", "GWR" (default), "MGWR" ,
#'   "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv", "MGWRSAR_0_kc_kv",
#'   "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details.
#' @param S_out  A matrix with variables used in kernel for out-sample prediction (jacknife mode)
#' @param mstop   Number of iterations for mboost.
#' @param nu  Learning rate for mboost.
#' @param noisland A boolean to avoid isle with no neighbours for non adaptive kernel, default FALSE
#' @return a list of objects for MGWRSAR wrapper
#' @noRd
GWR<-function(Y,XV,ALL_X,S,H,NN, kernels,adaptive=FALSE, Type = "GD",SE=FALSE, isgcv=FALSE,W=NULL,TP=NULL,dists=NULL,indexG=NULL,Wd=NULL,doMC=FALSE,ncore=1,Model,S_out=FALSE,mstop=150,nu=0.1,noisland=FALSE,get_ts=FALSE,family=NULL){
  SEV=NULL
  X=XV
  n<-nrow(ALL_X)
  if(is.null(Wd)){
    Z=S[TP,]
    if(S_out){
      pred=TRUE
      S<-S[-TP,]
      Y<-Y[-TP]
      XV<-XV[-TP,]
      ALL_X<-ALL_X[-TP,]
      NN=length(Y)
    } else pred=FALSE

      stage1=prep_w(H=H,kernels=kernels,coord_i=Z,coord_j=S,NN=NN,ncolX=ncol(XV),Type=Type,adaptive=adaptive,dists=dists,indexG=indexG,rowNorm=TRUE,noisland=noisland)
      indexG=stage1$indexG
      dists=stage1$dists
      Wd=stage1$Wd
    }
    if(Model %in% c('GWR_glmboost','GWR_gamboost_linearized')) model=gwr_beta_glmboost(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,doMC=doMC,ncore=ncore,pred=pred,mstop=mstop,nu=nu,family=family) else if (Model %in% c('GWR_glm'))  model=gwr_beta_glm(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,doMC=doMC,ncore=ncore,pred=pred,family=family) else model=gwr_beta(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,W=W,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,doMC=doMC,ncore=ncore,pred=pred,get_ts=get_ts)
  if(SE & !isgcv) list(Betav=model$Betav,SEV=model$SEV,edf=n-model$tS,tS=model$tS) else if(get_ts) list(Betav=model$Betav,SEV=NULL,edf=NULL,tS=model$tS)  else list(Betav=model$Betav,SEV=NULL,edf=NULL,tS=NULL)
}
