#' Estimation of linear and local linear model with spatial autocorrelation model (mgwrsar).
#'
#' MGWRSAR is is a wrapper function for estimating linear and local linear models
#' with spatial autocorrelation (SAR models with spatially varying coefficients).
#'
#' @usage MGWRSAR(formula,data,coords,fixed_vars=NULL,kernels,H,
#' Model='GWR',control=list())
#'
#' @param formula  a formula.
#' @param fixed_vars a vector with the names of spatiallay constant coefficient for
#' mixed model. All other variables present in formula are supposed to be spatially
#' varying. If empty or NULL (default), all variables in formula are supposed to be
#' spatially varying.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coords default NULL, a dataframe or a matrix with coordinates, not
#' required if data is a spatial dataframe.
#' @param Model character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' explanation.
#' @param kernels A vector containing the kernel types. Possible types:
#' rectangle ("rectangle"), bisquare ("bisq"), tricube ("tcub"), epanechnikov ("epane"), gaussian
#' ("gauss")) .
#' @param H vector containing the bandwidth parameters for the kernel functions.
#' @param control list of extra control arguments for MGWRSAR wrapper - see Details below
#' @details
#' \describe{
#' \item{Z}{ a matrix of variables for genralized kernel product, default NULL.}
#' \item{W}{ a row-standardized spatial weight matrix for Spatial Aurocorrelation, default NULL.}
#' \item{type}{ verbose mode, default FALSE.}
#' \item{adaptive}{A vector of boolean to choose adaptive version for each kernel.}
#' \item{kernel_w}{ the type of kernel for computing W, default NULL.}
#' \item{h_w}{ the bandwidth value for computing W, default 0.}
#' \item{Method}{ estimation technique for computing the models with Spatial Dependence. '2SLS' or 'B2SLS', default '2SLS'.}
#' \item{TP}{ A vector of target points, default NULL.}
#' \item{doMC}{Parallel computation, default FALSE}
#' \item{ncore}{number of CPU core for parallel computation, default 1}
#' \item{isgcv}{ computing LOOCV criteria (for example for selecting optimal bandwidth), default FALSE.}
#' \item{isfgcv}{ if TRUE, simplify the computation of CV criteria (remove or not i
#' when using local instruments for model with lambda spatially varying), default TRUE.}
#' \item{maxknn}{ when n >NmaxDist, only the maxknn first neighbours are used for distance compution, default 500.}
#' \item{NmaxDist}{ when n >NmaxDist only the maxknn first neighbours are used for distance compution, default 5000}
#' \item{verbose}{ verbose mode, default FALSE.}
#'}
#' @return MGWRSAR returns an object of class mgwrsar with at least the following components:
#' \describe{
#' \item{Betav}{ matrix of coefficients of dim(n,kv) x kv.}
#' \item{Betac}{ vector of coefficients of length kc.}
#' \item{Model}{ The sum of square residuals.}
#' \item{Y}{ The dependent variable.}
#' \item{XC}{ The explanatory variables with constant coefficients.}
#' \item{XV}{ The explanatory variables with varying coefficients.}
#' \item{X}{ The explanatory variables.}
#' \item{W}{ The spatial weight matrix for spatial dependence.}
#' \item{isgcv}{ if gcv has been computed.}
#' \item{edf}{ The estimated degrees of freedom.}
#' \item{formula}{The formula.}
#' \item{data}{ The dataframe used for computation.}
#' \item{Method}{ The type of model.}
#' \item{coords}{ The spatial coordinates of observations.}
#' \item{H}{ The bandwidth vector.}
#' \item{fixed_vars}{ The names of constant coefficients.}
#' \item{kernels}{ The kernel vector.}
#' \item{SSR}{ The sum of square residuals.}
#' \item{residuals}{ The vector of residuals.}
#' \item{fit}{ the vector of fitted values.}
#' \item{sev}{ local standard error of parameters.}
#' \item{get_ts}{Boolean, if trace of hat matrix Tr(S) should be stored.}
#' \item{NN}{ Maximum number of neighbors for weights computation}
#'}
#' MGWRSAR is is a wrapper function for estimating linear and local linear model
#' with spatial autocorrelation that  allows to estimate the following models :
#' \eqn{y=\beta_c X_c+\,\epsilon_i} (OLS)
#'
#' \eqn{y=\beta_v(u_i,v_i) X_v+\,\epsilon_i} (GWR)
#'
#' \eqn{y=\beta_c X_c+\beta_v(u_i,v_i) X_v+\,\epsilon_i} (MGWR)
#'
#' \eqn{y=\lambda Wy+\beta_c X_c+\,\epsilon_i} (MGWR-SAR(0,k,0))
#'
#' \eqn{y=\lambda Wy+\beta_v(u_i,v_i)X_v+\,\epsilon_i} (MGWR-SAR(0,0,k))
#'
#' \eqn{y=\lambda Wy+\beta_c X_c+\beta_v(u_i,v_i)X_v+\,\epsilon_i} (MGWR-SAR(0,k_c,k_v))
#'
#' \eqn{y=\lambda(u_i,v_i) Wy+\beta_c X_c+\,\epsilon_i} (MGWR-SAR(1,k,0))
#'
#' \eqn{y=\lambda(u_i,v_i)Wy+\beta_v(u_i,v_i)X_v+\,\epsilon_i} (MGWR-SAR(1,0,k))
#'
#' \eqn{y=\lambda(u_i,v_i)Wy+\beta_cX_c+\beta_v(u_i,v_i)X_v+\,\epsilon_i} (MGWR-SAR(1,k_c,k_v))
#'
#' When model imply spatial autocorrelation, a row normalized spatial weight matrix must be provided.
#'2SLS and Best 2SLS method can be used.
#' When model imply local regression, a bandwidth and a kernel type must be provided. Optimal bandwidth can be estimated
#' using bandwidths_mgwrsar function.
#' When model imply mixed local regression, the names of stationary covariates must be provided.
#'
#' #' In addition to the ability of considering spatial autocorrelation in GWR/MGWR like models,
#' MGWRSAR function introduces several useful technics for estimating local regression with space coordinates:
#' \itemize{
#' \item{it uses RCCP and RCCPeigen code that speed up computation and allows parallel computing via doMC package;}
#' \item{it allows to drop out variables with not enough local variance in local regression, which allows to consider dummies in GWR/MGWR framework without trouble.}
#' \item{it allows to drop out local outliers in local regression.}
#' \item{it allows to consider additional variable for kernel, including  time (asymetric kernel) and categorical variables (see Li and Racine 2010). Experimental version.}
#' }
#'
#' @references
#'
#' Geniaux, G. and Martinetti, D. (2017). A new method for dealing simultaneously with spatial autocorrelation and spatial heterogeneity in regression models. Regional Science and Urban Economics. (https://doi.org/10.1016/j.regsciurbeco.2017.04.001)
#'
#' McMillen, D. and Soppelsa, M. E. (2015). A conditionally parametric probit model of
#' microdata land use in chicago. Journal of Regional Science, 55(3):391-415.
#'
#' Loader, C. (1999). Local regression and likelihood, volume 47. springer New York.
#'
#' Franke, R. and Nielson, G. (1980). Smooth interpolation of large sets of scattered data.
#' International journal for numerical methods in engineering, 15(11):1691-1704.
#' @seealso  bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, predict_mgwrsar, kernel_matW
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix)
#'  ## of 4 nearest neighbors with 0 in diagonal
#'  W=kernel_matW(H=4,kernels='rectangle',coord_i=coords,NN=4,adaptive=TRUE,
#'  diagnull=TRUE,rowNorm=TRUE)
#'  mgwrsar_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,
#'  coords=coords, fixed_vars='X2',kernels=c('gauss'),H=20, Model = 'MGWRSAR_0_kc_kv',
#'  control=list(SE=FALSE,adaptive=TRUE,W=W))
#'  summary_mgwrsar(mgwrsar_0_kc_kv)
#' }
MGWRSAR <- function(formula, data, coords, fixed_vars = NULL, kernels, H,Model = "GWR", control = list()){
  set.seed(123)
 while(sum(duplicated(coords))>0) {
    coords<-jitter(coords,0.001)
    #warning('coords have been jittered because there is some duplicated location.')
 }
  colnames(coords)<-c('x','y')
    ptm<-proc.time()
    mycall <- match.call()
    n <-nrow(data)
    assign_control(control,n)
    gwrenv=environment()
    prep_var(gwrenv)
    if (Model == "OLS") {
      model<-list()
      lml <- lm.fit(as.matrix(X), as.matrix(Y))
      model$Betac <- lml$coefficients
      if (SE) {
        rss <- sum(lml$residuals^2)
        rdf <- length(Y) - ncol(X)
        resvar <- rss/rdf
        R <- chol2inv(lml$qr$qr)
        diagR=diag(R)
        model$se =  sqrt(diagR * resvar)
      }
      model$edf = rdf
      model$XC = X
      model$XV = NULL
    }
    else if (Model == "SAR") {
      model<-list()
      if (Method %in% c("B2SLS", "2SLS")) {
        keep = which(!is.na(coefficients(lm.fit(X, Y))))
        mymodel <- mod(as.matrix(Y), as.matrix(X[, keep]),
                       W, as.matrix(X), as.matrix(Y), rep(1, n), "L0",
                       (Method == "B2SLS"), FALSE, SE)
        Betac = mymodel$Betav
        if (SE)
          se = mymodel$se
        else se = NULL
        model$Betac[c(keep, ncol(X) + 1)] <- Betac
        if (SE) {
          model$se <- rep(0, ncol(X) + 1)
          model$se[c(keep, ncol(X) + 1)] <- se
        }
      }
      model$edf = n - ncol(X) - 1
      model$XC = as.matrix(cbind(X, as.matrix(W %*% Y)))
      model$XV = NULL
      model$Betav = NULL
      model$Y = Y
    }
    else if (Model %in%  c("GWR","GWR_glm","GWR_glmboost",'GWR_gamboost_linearized',"MGWRSAR_1_0_kv","GWR_multiscale") ){
      if(Model=="MGWRSAR_1_0_kv") Wx=W else Wx=NULL
      model = GWR(Y=Y,XV=XV,ALL_X=X,S=S,H=H,NN=NN,W=Wx, kernels=MykernelS,adaptive=adaptive, Type = Type,SE=SE, isgcv=isgcv,TP=TP,doMC=doMC,ncore=ncore,dists=dists,indexG=indexG,Wd=Wd,Model=Model,S_out=S_out,get_ts=get_ts,mstop=mstop,nu=nu,family=family)
      model$Betac = NULL
      model$XC = NULL
      if(Model=="MGWRSAR_1_0_kv") {
        if(is.null(new_data)) model$XV = cbind(XV, as.numeric(W %*% Y)) else model$XV = new_XV
        model$edf = model$edf - 1
      } else if(!S_out) model$XV = XV else model$XV = XV[TP,]
    } else {
      if(Model=="MGWR") Wx=NULL else Wx=W
      model<- MGWR(Y=Y,XC=XC,XV=XV,S=S,H=H,NN=NN, kernels=kernels,adaptive=adaptive, Type = Type,SE=SE, isgcv=isgcv,W=W,TP=TP,Model=Model,dists=dists,indexG=indexG,Wd=Wd,S_out=S_out,doMC=doMC,ncore=ncore,get_ts=get_ts)
      model$Y = Y
      model$XC=XC
      if(Model =="MGWRSAR_1_kc_kv") {
        if(is.null(new_data)) model$XV = cbind(XV, as.numeric(W %*% Y)) else {
          model$XV = new_XV
          model$XC = new_XC
        }
        model$edf = model$edf - 1
      } else if (Model =="MGWRSAR_1_kc_0") {
        if(is.null(new_data)) model$XV = as.numeric(W %*% Y) else {
          model$XV = NULL
          model$XC = new_XC
        }
      } else model$XV = XV
      if(Model=="MGWRSAR_0_kc_kv") {
        if(is.null(new_data)) model$XC = cbind(XC, as.numeric(W %*% Y)) else {
          model$XV = new_XV
          model$XC = new_XC
        }
        model$edf = model$edf - 1
      } else if(Model=="MGWRSAR_0_0_kv") {
        if(is.null(new_data)) model$XC = as.numeric(W %*% Y) else {
          model$XV = new_XV
          model$XC = NULL
        }
      }
    }
    term1 = 0
    term2 = 0
    if (!is.null(model$Betav) & is.null(new_data)) term1 <- rowSums(model$XV * model$Betav)
    if (!is.null(model$Betac) & is.null(new_data)) term2 <- model$XC %*% as.matrix(model$Betac)
    if(is.null(new_data)) {residuals <- Y - term1 - term2
    fit = as.numeric(term1 + term2)
    }

    if (!is.null(new_data)){
      if(!(Model %in% c('GWR','MGWR'))){
      if(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv'))  {
        lambda_pred=model$Betav[,ncol(model$Betav)]
        if(Model =='MGWRSAR_1_0_kv') beta_pred=model$Betav[,-ncol(model$Betav)] else beta_pred=cbind(matrix(model$Betac,nrow=nrow(model$Betav),ncol=length(model$Betac),byrow=TRUE),model$Betav[,-ncol(model$Betav)])
        } else if(Model %in% c('MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv')) {
          lambda_pred=model$Betac[length(model$Betac)]
          if(Model =='MGWRSAR_0_0_kv')  beta_pred=model$Betav else if(Model=='MGWRSAR_0_kc_kv') beta_pred=cbind(matrix(model$Betac[-length(model$Betac)],nrow=nrow(model$Betav),ncol=length(model$Betac)-1,byrow=TRUE),model$Betav)
        } else if(Model =='MGWRSAR_1_kc_0') {
          lambda_pred=model$Betav[,ncol(model$Betav)]
          beta_pred=matrix(model$Betac,nrow=nrow(model$Betav),ncol=length(model$Betac),byrow=TRUE)
        } else if(Model =='SAR') {
          lambda_pred=model$Betac[length(model$Betac)]
          beta_pred=model$Betac[-length(model$Betac)]
        }
        pred=list(lambda_pred=lambda_pred,beta_pred=beta_pred)
        } else  {
        term1 <- rowSums(new_XV * model$Betav)
        if(Model=='MGWR') term2 <- new_XC %*% as.matrix(model$Betac)
        pred = as.numeric(term1 + term2)
      }
    }

    try(colnames(model$Betav) <- names_betav, silent = TRUE)
    try(colnames(model$SEV) <- names_betav, silent = TRUE)
    try(names(model$Betac) <- names_betac, silent = TRUE)
    try(names(model$se) <- names_betac, silent = TRUE)
    if(!is.null(new_data)) {z <- list(Betav = model$Betav, Betac = model$Betac,pred = pred,XC = model$XC, XV = model$XV)} else {
      if(is.null(model$tS)) AICc=NULL else AICc<- n*log(sqrt(mean((residuals[TP])^2)))+n*log(2*pi)+n*(n+model$tS)/(n-2-model$tS)
      z <- list(Betav = model$Betav, Betac = model$Betac, sev = model$SEV,se = model$se, Model = Model, Y = Y, XC = model$XC, XV = model$XV,X = X, W = W, isgcv = isgcv, edf = model$edf, tS = model$tS,AICc=AICc,formula = formula, data = data, Method = Method, coords = coords,H = H, fixed_vars = fixed_vars, kernels = kernels,adaptive=adaptive, fit = fit,pred=pred,residuals = residuals, RMSE=sqrt(mean((residuals[TP])^2)),RMSEn=sqrt(mean((residuals)^2)),SSR = sum(residuals^2), Type = Type,S = S,NN=NN,TP=TP, doMC=doMC,ncore=ncore,mycall = mycall,ctime=(proc.time()-ptm)[3])
    class(z) <- "mgwrsar"}
    invisible(z)
  }
