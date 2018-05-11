#' Estimation of linear and local linear model with spatial autocorrelation model (mgwrsar).
#'
#' MGWRSAR is is a wrapper function for estimating linear and local linear models
#' with spatial autocorrelation (SAR models with spatially varying coefficients).
#'
#' @usage MGWRSAR(formula,data,coord,fixed_vars=NULL,kernels,H,
#' Model='GWR',control=list())
#'
#' @param formula  a formula.
#' @param fixed_vars a vector with the names of spatiallay constant coefficient for
#' mixed model. All other variables present in formula are supposed to be spatially
#' varying. If empty or NULL (default), all variables in formula are supposed to be
#' spatially varying.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coord default NULL, a dataframe or a matrix with coordinates, not
#' required if data is a spatial dataframe.
#' @param Model character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0". See Details for more
#' explanation.
#' @param kernels vector containing the kernel types. Possible types: k nearest
#' neighbors ("knn"), bisquare ("bisq"), adaptative bisquare ("bisq_knn"), gaussian
#' ("gauss"), adaptative gaussian ("gauss_adapt").
#' @param H vector containing the bandwidth parameters for the kernel funcitons.
#' @param control list of extra control arguments for MGWRSAR wrapper - see Details below
#' @details
#' \itemize{
#' \item{Z}{ a matrix of variables for genralized kernel product, default NULL.}
#' \item{W}{ a row-standardized spatial weight matrix for Spatial Aurocorrelation, default NULL.}
#' \item{type}{ verbose mode, default FALSE.}
#' \item{kernel_w}{ the type of kernel for computing W, default NULL.}
#' \item{h_w}{ the bandwidth value for computing W, default 0.}
#' \item{Method}{ estimation technique for computing the models with Spatial Dependence. '2SLS' or 'B2SLS', default '2SLS'.}
#' \item{isgcv}{ computing CV criteria (for example for selecting optimal bandwidth), default FALSE.}
#' \item{isfgcv}{ if TRUE, simplify the computation of CV criteria (remove or not i
#' when using local instruments for model with lambda spatially varying), default TRUE.}
#' \item{maxknn}{ when n >NmaxDist, only the maxknn first neighbours are used for distance compution, default 500.}
#' \item{NmaxDist}{ when n >NmaxDist only the maxknn first neighbours are used for distance compution, default 5000}
#' \item{verbose}{ verbose mode, default FALSE.}
#'}
#' @return MGWRSAR returns an object  of class mgwrsar with at least the following components:
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
#' \item{coord}{ The spatial coordinates of observations.}
#' \item{H}{ The bandwidth vector.}
#' \item{fixed_vars}{ The names of constant coefficients.}
#' \item{kernels}{ The kernel vector.}
#' \item{SSR}{ The sum of square residuals.}
#' \item{residuals}{ The vector of residuals.}
#' \item{fit}{ the vector of fitted values.}
#' \item{sev}{ local standard error of parameters.}
#'
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
#' @seealso  bandwidths_mgwrsar, summary_mgwrsar, plot_mgwrsar, predict_mgwrsar, kernelW_C
#' @examples
#' \donttest{
#' data(data_mgwrsar)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#' fixed_vars=NULL,kernels=c('gauss_knn'),
#' H=120, Model = 'GWR',control=list())
#' summary_mgwrsar(model_GWR)
#'
#' W=KNN(coord,8)
#' model_MGWRSAR_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,
#' coord=coord,fixed_vars='Intercept',kernels=c('gauss_adapt'),
#' H=120, Model = 'MGWRSAR_0_kc_kv',control=list(W=W))
#' summary_mgwrsar(model_MGWRSAR_0_kc_kv)
#' }
MGWRSAR <-
function(formula,data,coord,fixed_vars=NULL,kernels,H,Model='GWR',control=list())
{
mycall <- match.call()
n <- nrow(data)
########## control parameters
assign_control(control,n)
#### verifying parameters
if(is.null(coord)) { if(class(data) %in% c('SpatialPointsDataFrame','SpatialGridDataFrame','SpatialPixelsDataFrame')) coord=as.matrix(coordinates(data)) else stop("coord required")}

if(length(kernels)>1) SS=as.matrix(cbind(coord,Z)) else SS=as.matrix(coord)

if(is.null(H[1]) | is.null(kernels[1])) stop ("kernels list and bandwidths H required")

    if(is.null(fixed_vars) & Model %in% c('MGWR','MGWRSAR_0_kc_kv','MGWRSAR_1_kc_kv')) stop('You must provide fixed_vars for mixed models')
    if(is.null(W) & Model %in% c('SAR','MGWRSAR_1_0_kv','MGWRSAR_0_0_kv','MGWRSAR_0_kc_kv','MGWRSAR_1_kc_kv','MGWRSAR_1_kc_0')) stop('You must provide W for models with spatial dependence')
    if(!is.null(fixed_vars) & Model %in% c('GWR','SAR','MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')) {
    	fixed_vars=NULL
    	if(verbose) cat('\n-----------------------------------------------------\nfixed_vars set to NULL because model= ',Model,'\n-----------------------------------------------------\n')
    	}
if(!is.null(W) & Model %in% c('GWR','OLS','MGWR')) {
    	if(verbose) cat('\n-----------------------------------------------------\nW not used because model= ',Model,'\n-----------------------------------------------------\n')
    	}

    ####
	mf <- model.frame(formula,data)
	mt <- attr(x = mf, which = "terms")
	X=model.matrix(object = mt, data = mf, contrasts.arg = contrasts)
    Y <- model.extract(mf, "response")
    idx1 <- match("(Intercept)", colnames(X))
    if (!is.na(idx1))
        colnames(X)[idx1] <- "Intercept"
    if (!is.null(fixed_vars)) {
        idx.fixed <- match(fixed_vars, colnames(X))
        XC <- as.matrix(X[, idx.fixed])
        colnames(XC) <- colnames(X)[idx.fixed]
        if (length(idx.fixed) < ncol(X)) {
            XV <- as.matrix(X[, -idx.fixed])
            colnames(XV) <- colnames(X)[-idx.fixed]
        }
        else XV = NULL
    }
    else {
        XV = as.matrix(X)
        XC = NULL
    }
    Y <- as.matrix(Y)
    model <- c()
    coord = as.matrix(coord)
    ## W vide pour MGWR
    if(is.null(W)) W<-as(Matrix(0,nrow=n,ncol=n),'dgCMatrix')
    #### nom variables
    names_betac= colnames(XC)
    names_betav= colnames(XV)
    if(Model %in% c('OLS')) names_betac=colnames(X)
    if(Model %in% c('SAR')) names_betac=c(colnames(X),'lambda')
    if(Model %in% c('MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv')) names_betac=c(names_betac,'lambda')
    if(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv')) names_betav=c(names_betav,'lambda')
    if(Model=='MGWRSAR_1_kc_0') {names_betav=c('lambda');names_betac=colnames(X);}

    ## C variables
    MykernelS=kernels
	HH=H
	Y=as.matrix(Y)
	X=as.matrix(X)
	if(!is.null(XC)) XC=as.matrix(XC)
	if(!is.null(XV)) XV=as.matrix(XV)

	##if(Model %in% c('MGWRSAR_1_0_kv','MGWRSAR_0_0_kv')) XC=rep(1,n)
    ###
    if (Model == "OLS") {
        mymodel<-fastlmLLT_C(as.matrix(X),as.matrix(Y),SE)
        model$Betac <- mymodel$Betav
        if(SE)  model$se= mymodel$se
        model$edf=n-ncol(X)
        model$XC=X
        model$XV=NULL
    } else if (Model == "SAR") {
        if (Method %in% c("B2SLS", "2SLS")) {
        keep=which(!is.na(coefficients(lm.fit(X,Y))))
        mymodel<- mod(as.matrix(Y),as.matrix(X[,keep]),W,as.matrix(X),as.matrix(Y),rep(1,n),'L0',(Method =="B2SLS"), FALSE,SE)
        Betac = mymodel$Betav
        if(SE) se= mymodel$se else se=NULL
        model$Betac<-rep(0,ncol(X)+1)
        model$Betac[c(keep,ncol(X)+1)]<-Betac
        if(SE) {
            model$se<-rep(0,ncol(X)+1)
            model$se[c(keep,ncol(X)+1)]<-se
          }
        }
        model$edf=n-ncol(X)-1
        model$XC=as.matrix(cbind(X,as.matrix(W%*%Y)))
        model$XV=NULL
        model$Y=Y
    } else if (Model == "GWR"){
     model=GWR(Y,XV,X,Y,S=SS,H=H,kernels=MykernelS,type = Type, minv = minv, maxknn = maxknn, NmaxDist = NmaxDist, SE=SE,isgcv=isgcv, TIME=TIME, decay=decay,W=NULL,betacor=Lambdacor,remove_local_outlier=remove_local_outlier,outv=outv,doMC=doMC,ncore=ncore,Wh=Wh,xratiomin=xratiomin)
     model$Betac=NULL
     model$XV=XV
     model$XC=NULL
    } else if (Model == "MGWRSAR_1_0_kv"){
    model=GWR(Y,XV,X,Y,S=SS,H=H,kernels=MykernelS,type = Type, minv = minv, maxknn = maxknn, NmaxDist = NmaxDist, SE=SE, isgcv=isgcv, TIME=TIME, decay=decay,W=W,betacor=Lambdacor,remove_local_outlier=remove_local_outlier,outv=outv,doMC=doMC,ncore=ncore,Wh=Wh,xratiomin=xratiomin)
     model$Betac=NULL
     model$XV=cbind(XV,as.numeric(W%*%Y))
     model$XC=NULL
     model$edf=model$edf-1
     } else {
    model <- MGWR(Y,XC,XV,S=SS,H=H,kernels=MykernelS,type = Type,model=Model, minv = minv, maxknn = maxknn, NmaxDist = NmaxDist, SE=SE, isgcv=isgcv, TIME=TIME, decay=decay,W=W,betacor=Lambdacor,remove_local_outlier=remove_local_outlier,outv=outv,doMC=doMC,ncore=ncore,Wh=Wh,xratiomin=10e-10)
    model$Y=Y
    }


    ####### Fin estim
    term1 = 0
    term2 = 0
    if (!is.null(model$Betav))
        term1 <- rowSums(model$XV * model$Betav)
    if (!is.null(model$Betac))
        term2 <- model$XC %*% as.matrix(model$Betac)
    residuals <- Y - term1 - term2
    fit=Y-residuals
    try(colnames(model$Betav)<-names_betav,silent =TRUE)
    try(colnames(model$SEV)<-names_betav,silent =TRUE)
    try(names(model$Betac)<-names_betac,silent =TRUE)
    try(names(model$se)<-names_betac,silent =TRUE)
    z <- list(Betav = model$Betav, Betac = model$Betac, sev=model$SEV, se=model$se, Model = Model,
        Y = Y, XC = model$XC, XV =  model$XV,X=X, W = W, isgcv = isgcv,edf=model$edf,tS=model$tS,formula=formula,data=data,Method=Method,coord=coord,H=H,fixed_vars=fixed_vars,kernels=kernels,fit=fit,residuals=residuals,SSR=sum(residuals^2),type = Type,S=SS,mycall=mycall)
    class(z)<-'mgwrsar'
    #rm(envavirer,envir=.GlobalEnv)
    invisible(z)
}
