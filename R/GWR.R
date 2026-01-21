#' Internal GWR and MGWRSAR Estimation Wrapper
#'
#' @description
#' This function acts as the core backend for estimating Geographically Weighted Regression (GWR)
#' and various Mixed GWR-SAR models (MGWRSAR). It handles the preparation of spatial weight matrices
#' (stage 1) and dispatches the estimation to specific sub-routines based on the selected `Model`
#' and `family`.
#'
#' @references
#' Geniaux, G. and Martinetti, D. (2018). A new method for dealing with spatial heterogeneity in coefficients of spatial models. Regional Science and Urban Economics.
#' Geniaux, G. (2024). Speeding up MGWRSAR.
#'
#' @param Y A numeric vector of the response variable.
#' @param XV A matrix or data frame with covariates assumed to have stationary parameters (Global part).
#' @param ALL_X A matrix or data frame containing all covariates (both spatially varying `XC` and stationary `XV`).
#' @param S A matrix containing the spatial coordinates (or variables used in the kernel).
#' @param H A vector of bandwidths.
#' @param NN An integer indicating the number of spatial neighbours for kernel computations.
#' @param kernels A character vector specifying the kernel types (e.g., 'gauss', 'bisq').
#' @param adaptive A logical vector or boolean; if TRUE, adaptive bandwidths (nearest neighbours) are used. Default is FALSE.
#' @param Type A character string specifying the type of Generalized Kernel product.
#'   Options include: 'GD' (spatial only), 'GDC' (spatial + categorical), 'GDX' (spatial + continuous), etc.
#' @param SE Logical; if TRUE, standard errors are computed. Default is FALSE.
#' @param isgcv Logical; if TRUE, leave-one-out cross-validation is performed. Default is FALSE.
#' @param W A spatial weight matrix (sparse or dense) for spatial autocorrelation (SAR component). Default is NULL.
#' @param TP A vector of indices for target points (where the estimation is performed). Default is NULL.
#' @param dists A precomputed matrix of spatial distances. Default is NULL.
#' @param indexG A precomputed matrix of indices of Nearest Neighbours. Default is NULL.
#' @param Wd A precomputed matrix of spatial weights. Default is NULL.
#' @param ncore Integer; the number of cores for parellization.
#' @param Model A character string specifying the type of model to estimate.
#'   Possible values: "OLS", "SAR", "GWR" (default), "MGWR",
#'   "MGWRSAR_0_0_kv", "MGWRSAR_1_0_kv", "MGWRSAR_0_kc_kv",
#'   "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0".
#' @param TP_estim_as_extrapol Logical; if TRUE, treats target points as extrapolation points (excludes them from calibration).
#' @param mstop Integer; the number of boosting iterations (for `GWR_glmboost` models). Default is 150.
#' @param nu Numeric; the learning rate for boosting. Default is 0.1.
#' @param noisland Logical; if TRUE, prevents formation of islands with no neighbours for non-adaptive kernels. Default is FALSE.
#' @param get_ts Logical; if TRUE, returns the Hat matrix (Trace S). Default is FALSE.
#' @param get_s Logical; if TRUE, returns the full S matrix. Default is FALSE.
#' @param get_Rk Logical; if TRUE, returns the resolution matrix Rk. Default is FALSE.
#' @param family An object of class `family` (e.g., `gaussian()`, `binomial()`). Used for GLM-based GWR.
#' @param alpha Numeric; a scaling or adjustment parameter for the weights (depending on the kernel type). Default is 1.
#'
#' @return A list containing the following components (depending on the input parameters):
#'   \item{Betav}{Matrix of estimated coefficients.}
#'   \item{SEV}{Matrix of standard errors (if `SE=TRUE`).}
#'   \item{edf}{Effective degrees of freedom.}
#'   \item{tS}{Trace of the Hat matrix.}
#'   \item{Shat}{The Hat matrix (if requested).}
#'   \item{Rk}{Resolution matrix (if requested).}
#'
#' @keywords internal
#' @noRd
GWR<-function(Y,XV,ALL_X,S,H,NN, kernels,adaptive=FALSE, Type = "GD",SE=FALSE, isgcv=FALSE,W=NULL,TP=NULL,dists=NULL,indexG=NULL,Wd=NULL,ncore=1,Model,TP_estim_as_extrapol=FALSE,mstop=150,nu=0.1,noisland=FALSE,get_ts=FALSE,get_s=FALSE,get_Rk=FALSE,family=NULL,alpha=1){
  SEV=NULL
  X=XV
  n<-nrow(ALL_X)
  if(is.null(Wd)){
    Z=S[TP,]
    if(TP_estim_as_extrapol){
     Y<-Y[-TP]
     XV<-XV[-TP,]
     ALL_X<-ALL_X[-TP,]
    }
    stage1=prep_w(H=H,kernels=kernels,Type=Type,adaptive=adaptive,dists=dists,indexG=indexG,alpha=alpha)
      indexG=stage1$indexG
      dists=stage1$dists
      Wd=stage1$Wd
      isolated_idx <- integer(0)

      if (kernels[1] != "gauss") {
        isolated_idx <- which(rowSums(Wd > 0) < ncol(XV))
        if(length(isolated_idx)>0){
            assign('isolated_idx',isolated_idx,envir=parent.frame())
        }
      }

  }
    if(Model %in% c('GWR_glmboost','GWR_gamboost_linearized')) model=gwr_beta_glmboost(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,ncore=ncore,TP_estim_as_extrapol=TP_estim_as_extrapol,mstop=mstop,nu=nu,family=family) else if (family$family %in% c("binomial","quasibinomial"))  model=gwr_beta_glm(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,ncore=ncore,TP_estim_as_extrapol=TP_estim_as_extrapol,family=family,get_ts=get_ts,get_s=get_s) else {
      if(length(TP)==n) model=gwr_beta_pivotal_qrp_full(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,W=W,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,ncore=ncore,TP_estim_as_extrapol=TP_estim_as_extrapol,get_ts=get_ts,get_s=get_s,get_Rk=get_Rk,isolated_idx=NULL) else model=gwr_beta(Y=Y,XV=XV,ALL_X=ALL_X,TP=TP,indexG=indexG,Wd=Wd,NN=NN,W=W,isgcv=isgcv,SE=SE,H=H,kernels=kernels,adaptive=adaptive,ncore=ncore,TP_estim_as_extrapol=TP_estim_as_extrapol,get_ts=get_ts,get_s=get_s,get_Rk=get_Rk,isolated_idx=NULL)
}

  if(SE & !isgcv) list(Betav=model$Betav,SEV=model$SEV,edf=n-model$tS,tS=model$tS,Shat=model$Shat,TS=model$TS,Rk=model$Rk) else if(get_ts) list(Betav=model$Betav,SEV=NULL,edf=NULL,tS=model$tS,Shat=model$Shat,TS=model$TS,Rk=model$Rk)  else list(Betav=model$Betav,SEV=NULL,edf=NULL,tS=NULL,Shat=NULL,TS=NULL)
}
