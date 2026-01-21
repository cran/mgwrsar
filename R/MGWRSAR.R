#' Estimation of linear and local linear model with spatial autocorrelation model (mgwrsar).
#'
#' MGWRSAR is is a wrapper function for estimating linear and local linear models
#' with spatial autocorrelation (SAR models with spatially varying coefficients).
#'
#' @usage MGWRSAR(formula, data, coords, fixed_vars = NULL, kernels, H,
#' Model = "GWR", control = list())
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
#' \item{Z}{A matrix of variables for genralized kernel product, default NULL.}
#' \item{W}{A row-standardized spatial weight matrix for Spatial
#'  Aurocorrelation, default NULL.}
#' \item{Type}{Verbose mode, default FALSE.}
#' \item{adaptive}{A vector of boolean to choose adaptive version for
#'  each kernel.}
#' \item{kernel_w}{The type of kernel for computing W, default NULL.}
#' \item{h_w}{The bandwidth value for computing W, default 0.}
#' \item{Method}{Estimation method for computing the models with Spatial
#' Dependence. '2SLS' or 'B2SLS', default '2SLS'.}
#' \item{TP}{Avector of target points, default NULL.}
#' \item{ncore}{Number of CPU core for parallel computation, default 1}
#' \item{isgcv}{If TRUE, compute a LOOCV criteria, default FALSE.}
#' \item{maxknn}{When n >NmaxDist, only the maxknn first neighbours are used
#' for distance compution, default 500.}
#' \item{NmaxDist}{When n >NmaxDist only the maxknn first neighbours are used
#' for distance compution, default 5000}
#' \item{verbose}{Verbose mode, default FALSE.}
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
#' \item{fixed_vars}{The names of constant coefficients.}
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
#' using goldens_search_bandwid function.
#' When model imply mixed local regression, the names of stationary covariates must be provided.
#'
#' #' In addition to the ability of considering spatial autocorrelation in GWR/MGWR like models,
#' MGWRSAR function introduces several useful technics for estimating local regression with space coordinates:
#' \itemize{
#' \item{it uses RCCP and RCCPeigen code that speed up computation and allows parallel computings;}
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
#' @seealso  golden_search_bandwidth, summary, plot, predict, kernel_matW
#' @examples
#' \donttest{
#'  library(mgwrsar)
#'  ## loading data example
#'  data(mydata)
#'  coords=as.matrix(mydata[,c("x","y")])
#'  ## Creating a spatial weight matrix (sparce dgCMatrix)
#'  ## of 4 nearest neighbors with 0 in diagonal
#'  W=kernel_matW(H=4,kernels='rectangle',coords=coords,NN=4,adaptive=TRUE,
#'  diagnull=TRUE)
#'  mgwrsar_0_kc_kv<-MGWRSAR(formula = 'Y_mgwrsar_0_kc_kv~X1+X2+X3', data = mydata,
#'  coords=coords, fixed_vars='X2',kernels=c('gauss'),H=20, Model = 'MGWRSAR_0_kc_kv',
#'  control=list(SE=FALSE,adaptive=TRUE,W=W))
#'  summary(mgwrsar_0_kc_kv)
#' }
MGWRSAR <- function(formula, data, coords, fixed_vars = NULL, kernels, H,Model = "GWR", control = list()){
  set.seed(123, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion")

  control <- .mgwrsar_normalize_parallel_control(control, context = "MGWRSAR")
  get_SFdata()

  mf <- model.frame(formula, data)
  data<-data[,names(mf)]
  if(is.null(control$Type)) control$Type='GD'
  if(Model!='OLS'){
  if(control$Type %in% c('GD','GDT')) coords<-make_unique_by_structure(coords)
  if(control$Type %in% c('GDT','T')) control$Z<-make_unique_by_structure(control$Z)
  }

  # if(control$Type=='GDT'){
  #   if(kernels[1]!='gauss' | unlist(str_split(kernels[2], '_'))[1]!='gauss' ) stop('For Type=GDT only gauss kernel should be used')
  # }
    #colnames(coords)<-c('x','y')
  if (is.null(coords)) stop("coords must be provided")
  if (is.null(control$Type)) control$Type='GD'

    type_local <- if (!is.null(control$Type)) control$Type else "GD"

    if (type_local == "T") {
      if (is.null(colnames(coords))) colnames(coords) <- "t"
    } else {
      if (ncol(coords) == 2 && is.null(colnames(coords))) {
        colnames(coords) <- c("x", "y")
      } else if (ncol(coords) > 2) {
        colnames(coords)[1:2] <- c("x", "y")
      }
    }
    ptm<-proc.time()
    formula<-formula(lm(formula,data))
    mycall <- match.call()
    n <-nrow(data)
    assign_control(control,n)
    if(is.null(control$family)){
      control$family<-family<-gaussian(link = "identity")
    }
    gwrenv=environment()
    prep_var(gwrenv)
    if (Model == "OLS") {
      model<-list()
      lml <- lm.fit(as.matrix(X), as.matrix(Y))
      model$Betac <- lml$coefficients
      if(sum(is.na(model$Betac))>0) stop(paste0(names(model$Betac)[is.na(model$Betac)],' are fully collinear, remove these terms from formula'))
      if (SE) {
        rss <- sum(lml$residuals^2)
        rdf <- length(Y) - ncol(X)
        resvar <- rss/rdf
        R <- chol2inv(lml$qr$qr)
        diagR=diag(R)
        model$se =  sqrt(diagR * resvar)
        model$edf = rdf
      }

      model$XC = X
      model$tS=ncol(X)
      model$TS=lm.influence(lml)$hat
      if(get_s) {
        XXtX<-solve(crossprod(X)) %*% t(X)
        model$Shat=X%*% XXtX
        }
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
          names(model$se)<-names_betac[-length(names_betac)]
        }
      }
      model$edf = n - ncol(X) - 1
      model$XC = as.matrix(cbind(X, as.matrix(W %*% Y)))
      model$XV = NULL
      model$Betav = NULL
      model$Y = Y
      names(model$Betac)<-names_betac
    }
    else if (Model %in%  c("GWR","multiscale_gwr","GWR_glm","GWR_glmboost",'GWR_gamboost_linearized',"MGWRSAR_1_0_kv","GWR_multiscale") ){
      if(Model=="MGWRSAR_1_0_kv") Wx=W else Wx=NULL
      model = GWR(Y=Y,XV=XV,ALL_X=X,S=S,H=H,NN=NN,W=Wx, kernels=MykernelS,adaptive=adaptive, Type = Type,SE=SE, isgcv=isgcv,TP=TP,ncore=ncore,dists=dists,indexG=indexG,Wd=Wd,Model=Model,TP_estim_as_extrapol=TP_estim_as_extrapol,get_ts=get_ts,get_s=get_s,get_Rk=get_Rk,mstop=mstop,nu=nu,family=family,alpha=alpha)

      ## MODEL cases
      model$Betac = NULL
      model$XC = NULL
      model$XV=XV
      if(Model=="MGWRSAR_1_0_kv") {
        if(is.null(new_data)) model$XV = cbind(XV, as.numeric(W %*% Y))
        model$edf = model$edf - 1
      }  #else if(!is.null(TP_estim_as_extrapol)) model$XV = XV else model$XV = XV[TP,]
    } else {
      if(Model=="MGWR") Wx=NULL else Wx=W
      model<- MGWR(Y=Y,XC=XC,XV=XV,S=S,H=H,NN=NN, kernels=kernels,adaptive=adaptive, Type = Type,SE=SE, isgcv=isgcv,W=W,TP=TP,Model=Model,dists=dists,indexG=indexG,Wd=Wd,TP_estim_as_extrapol=TP_estim_as_extrapol,ncore=ncore,get_ts=get_ts,get_s=get_s,get_Rk=get_Rk,alpha=alpha)
      model$Y = Y
      if(!is.null(XC)) model$XC=XC
      if(Model =="MGWRSAR_1_kc_kv") {
        if(is.null(new_data)) model$XV = cbind(XV, as.matrix(W %*% Y,ncol=1))
        model$edf = model$edf - 1
      } else if (Model =="MGWRSAR_1_kc_0") {
        if(is.null(new_data)) model$XV = as.matrix(W %*% Y,ncol=1)
      } else model$XV = XV
      if(Model=="MGWRSAR_0_kc_kv") {
        if(is.null(new_data)) model$XC = cbind(XC, as.matrix(W %*% Y,ncol=1))
        model$edf = model$edf - 1
      } else if(Model=="MGWRSAR_0_0_kv") {
        if(is.null(new_data)) model$XC = as.matrix(W %*% Y,ncol=1)
      }
    }

    if(length(TP)<length(Y) & !TP_estim_as_extrapol & !sum(is.na(model$Betav)>0)){
      ## version Wtp
      if(!(kernels[1] %in% c('epane','bisq','triangle','tcub'))) {
      #Wtp=kernel_matW(H=H,kernels=kernels,coords=S[c(TP,(1:n)[-TP]),],NN=NN,TP=(1:length(TP)),Type=Type,adaptive=adaptive,diagnull=FALSE,alpha=alpha,dists=NULL,indexG=NULL,extrapol=F,QP=NULL)[,-(1:length(TP))]
      Wtp=kernel_matW(H=H,kernels=kernels,coords=S,NN=NN,TP=TP,Type=Type,adaptive=adaptive,diagnull=FALSE,alpha=alpha,dists=NULL,indexG=NULL,extrapol=F,QP=NULL)[,-TP]
      Wtp<- normW(Matrix::t(Wtp))
      } else {
      ## version shepard
      Wtp=kernel_matW(H=8,kernels='shepard',coords=S,NN=8+2,TP=(1:n)[-TP],Type=Type,adaptive=FALSE,diagnull=FALSE,extrapol=T,QP=TP)
      }

      model$Betav[-TP,]=as.matrix(Wtp%*% model$Betav[TP,])
      TPTS_approx<-lm.influence(lm.fit(X,Y))$hat
      model$TS[-TP]<-TPTS_approx[-TP]
      model$tS<-sum(model$TS)
      if(SE) model$SEV[-TP,]=as.matrix(Wtp%*% model$SEV[TP,])
      if(get_s){
        S<-matrix(0,nrow=n,ncol=n)
        S[TP,]<-model$Shat
        Z=as.matrix(Wtp%*% model$Shat)
        for(k in 1:ncol(XV)) S[-TP,]<-S[-TP,]+ Z*XV[-TP,k]
        model$Shat=S
      }

    }

    ### output
    mymodel<-new('mgwrsar')
    if(!is.null(model$TS)) mymodel@TS=as.numeric(model$TS)
    mymodel@Y = as.numeric(Y)
    mymodel@X = X
    if(!is.null(model$XC)) mymodel@XC = model$XC
    if(!is.null(model$XV)) mymodel@XV = model$XV
    mymodel@Model = Model
    if(!is.null(model$XC)) mymodel@XC = as.matrix(model$XC)
    if(!is.null(model$XV)) mymodel@XV =  as.matrix(model$XV)
    if(!is.null(W)) if(!is.null(dim(W))) mymodel@W = W

    mymodel@formula = as.formula(formula)
    mymodel@data = data[,names(mf)]
    mymodel@Method = Method
    mymodel@coords = coords
    if(!is.null(control$Z)) mymodel@Z=control$Z
    if(!is.null(fixed_vars)) mymodel@fixed_vars=fixed_vars
    if(!is.null(kernels)) mymodel@kernels=kernels
    mymodel@adaptive=adaptive
    mymodel@Type=Type
    if(Type=='GDT') {
      mymodel@alpha=alpha
      mymodel@Ht=H[2]
    }
    if(Type=='T') {
      mymodel@alpha=alpha
      mymodel@Ht=H[1]
    } else  mymodel@H=H[1]

    if(!is.null(NN)) mymodel@NN = NN
    mymodel@TP=TP
    mymodel@ncore=ncore
    mymodel@mycall=mycall
    mymodel@ctime=(proc.time()-ptm)[3]

    gwrenv=environment()
    mymodel<-format_and_diagno(e=gwrenv)

    invisible(mymodel)
  }
