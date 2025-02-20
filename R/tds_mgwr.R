#' Top-Down Scaling approach of multiscale GWR
#'
#' This function performs a multiscale Geographically Weighted Regression
#' (GWR) using a top-down scaling approach, adjusting GWR coefficients with
#' a progressively decreasing bandwidth as long as the AICc criterion improves.
#'
#' @usage tds_mgwr(formula,data,coords,Model='tds_mgwr',kernels='triangle',
#' fixed_vars=NULL,H2=NULL,control_tds=list(nns=30,get_AIC=FALSE),
#' control=list(adaptive=TRUE))
#' @param formula  a formula.
#' @param data a dataframe.
#' @param coords default NULL, a dataframe or a matrix with coordinates
#' @param Model character containing the type of model: Possible values are
#' "tds_mgwr" and "atds_mgwr", See Details for more explanation.
#' @param kernels A vector containing the kernel types. Possible types:
#' triangle ("triangle"), rectangle ("rectangle"), bisquare ("bisq"),
#' tricube ("tcub"), gaussian ("gauss"), epanechnikov ("epane").
#' @param H2 A scalar or vector of time bandwidths.
#' @param fixed_vars a vector with the names of spatiallay constant
#' coefficient for mixed model. All other variables present in formula
#' are supposed to be spatially varying. If empty or NULL (default),
#' all variables in formula are supposed to be spatially varying.
#' @param control_tds list of extra control arguments for tds_mgwr models
#' @details
#' \describe{
#' \item{nns}{Length of the sequence of decreasing bandwidth. Should be
#' between 20 and 100, default 30}
#' \item{get_AIC}{Boolean, if the Global AICc using Yu et al 2019 should be
#' computed. Required if the second stage 'atds_mgwr' has to be estimated.
#' default FALSE}
#' \item{init_model}{Starting model, 'GWR' or 'OLS', 'default OLS'.}
#' \item{model_stage1}{If model='tds_mgwr', model_stage1 can be used as a
#' starting model (either a GWR model or a preious tds_mgwr model).
#' For model='atds_mgwr, the user can specified an tds_mgwr
#' model already computed with get_AIC=TRUE. default NULL.}
#' \item{doMC}{Parallel computation, default FALSE.}
#' \item{ncore}{number of CPU core for parallel computation, default 1}
#' \item{tol}{Tolerance for stopping criteria, default 0.0001}
#' \item{nrounds}{Number of nrounds for 'atds_mgwr' model. Default 3.}
#' \item{verbose}{verbose mode, default FALSE.}
#' \item{V}{A vector of decreasing bandwidths given by the user, default NULL}
#' \item{first_nn}{The value of the highest bandwidth for the sequence of
#' decreasing bandwidth, default NULL.}
#' \item{minv}{The value of the smallest bandwidth for the sequence of
#' decreasing bandwidth, default number of covariates + 2 . }
#' \item{H}{A vector of bandwidth, default NULL}
#' }
#' @param control list of extra control arguments for MGWRSAR wrapper
#' @details
#' \describe{
#' \item{Z}{A matrix of variables for genralized kernel product, default NULL.}
#' \item{W}{A row-standardized spatial weight matrix for Spatial
#'  Aurocorrelation, default NULL.}
#' \item{type}{Verbose mode, default FALSE.}
#' \item{adaptive}{A vector of boolean to choose adaptive version for
#'  each kernel.}
#' \item{kernel_w}{The type of kernel for computing W, default NULL.}
#' \item{h_w}{The bandwidth value for computing W, default 0.}
#' \item{Method}{Estimation method for computing the models with Spatial
#' Dependence. '2SLS' or 'B2SLS', default '2SLS'.}
#' \item{TP}{Avector of target points, default NULL.}
#' \item{doMC}{Parallel computation, default FALSE. If TRUE and
#'  control_tds$doMC is also TRUE, then control$doMC is set to FALSE.}
#' \item{ncore}{Number of CPU core for parallel computation, default 1}
#' \item{isgcv}{If TRUE, compute a LOOCV criteria, default FALSE.}
#' \item{isfgcv}{If TRUE, simplify the computation of CV criteria
#' (remove or not i when using local instruments for model with lambda
#' spatially varying), default TRUE.}
#' \item{maxknn}{When n >NmaxDist, only the maxknn first neighbours are used
#' for distance compution, default 500.}
#' \item{NmaxDist}{When n >NmaxDist only the maxknn first neighbours are used
#' for distance compution, default 5000}
#' \item{verbose}{Verbose mode, default FALSE.}
#' }
#' @seealso  gwr_multiscale, MGWRSAR, bandwidths_mgwrsar, summary_mgwrsar.
tds_mgwr<-function(formula,data,coords,Model='tds_mgwr',kernels='triangle',fixed_vars=NULL,H2=NULL,control_tds=list(nns=30,get_AIC=FALSE),control=list(adaptive=TRUE)){
  start<-proc.time()
  init_param_tds()
  if(is.null(control_tds$model_stage1)){
    built_Vseq()
  } else if( model_stage1@Model %in% c('OLS','GWR')) {
    built_Vseq()
    control_tds$V<-control_tds$model_stage1@V<-model_stage1@V<-V
  }
  if(is.null(control_tds$H)){
         if(is.null(control_tds$model_stage1)){
        if(!is.null(TRUEBETA)) init_RMSE_history()
        stage1_tds_mgwr()
      }
      if(Model=='atds_mgwr'){
        stage2_atds_mgwr()
      }
    } else {
      if(!is.null(TRUEBETA)) init_RMSE_history()
      stage1_tds_mgwr_H()
    }
  returned_model
}

