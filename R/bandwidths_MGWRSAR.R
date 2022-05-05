#' bandwidths_mgwrsar
#' Select optimal kernel and bandwidth from a list of models, kernels and bandwidth candidates.
#'
#' Given a lm formula and a dataframe with coordinates, function bandwidths_mgwrsar optimizes the choice of
#' a bandwidth value for each of the chosen models and kernel types using a leave-one-out cross validation criteria.
#' A cross validated criteria is also used for selecting the best kernel type for a given model.
#'
#' @usage bandwidths_mgwrsar(formula, data,coord,
#' fixed_vars='Intercept',Models='GWR',Kernels_candidates='bisq',
#' control=list(),control_search=list())
#'
#' @param formula  a formula.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coord a dataframe or a matrix with coordinates, not required if data is a spatial dataframe, default NULL.
#' @param fixed_vars a vector with the names of spatially constant coefficient. For mixed model, if NULL, the default
#' #' is set to 'Intercept'.
#' @param Models character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0".
#' @param Kernels_candidates a vector with the names of kernel type.
#' @param control list of extra control arguments for MGWRSAR wrapper - see MGWRSAR help.
#' @param control_search list of extra control arguments for bandwidth/kernel search - see section below.
#'  @details
#'  \itemize{
#' \item{search_W}{if TRUE select an optimal spatial weight matrix using a moment estimator, default FALSE.}
#' \item{kernels_w}{if search_W is TRUE, kernels_w is a vector of candidated kernels types, default NULL.}
#' \item{lower_c}{lower bound for bandwidth search (default, the approximate first decile of distances).}
#' \item{upper_c}{upper bound for bandwidth search  (default, the approximate last decile of distances).}
#' \item{lower_d}{lower bound for discrete kernels, default 2*k+1.}
#' \item{lower_dW}{ower bound for discrete kernels for finding optimal spatial weight matrix, default 2.}
#' \item{lower_cW}{lower bound for  bandwidth search for finding optimal spatial
#' weight matrix  (default approximate 0.005 quantile of distances).}
#'}
#' @return bandwiths_MGWRSAR returns  a list with:
#' \item{config_model}{a vector with information about model, optimal kernel and
#' bandwidth for local regression, and optimal kernel and bandwith for spatial weight matrix W.}
#' \item{SSR}{The sum of square residuals.}
#' \item{CV}{The CV criteria.}
#' \item{model}{objects of class mgwrsar estimated using config_model}
#'
#' @details  Given a lm formula and a dataframe with coordinates, for each model in
#'  \code{Models} for wich a bandwidth is required, this function optimizes the choice of
#' a bandwidth value for each of the chosen models and kernel types using a leave one out
#' cross validation criteria. A cross validated criteria is also used for selecting the best kernel type for a given model.
#'
#' @references
#'
#'Geniaux, G. and Martinetti, D. (2017). A new method for dealing simultaneously with spatial autocorrelation and spatial heterogeneity in regression models. Regional Science and Urban Economics. (https://doi.org/10.1016/j.regsciurbeco.2017.04.001)
#'
#'
#' McMillen, D. and Soppelsa, M. E. (2015). A conditionally parametric probit model of
#' microdata land use in chicago. Journal of Regional Science, 55(3):391-415.
#'
#' Loader, C. (1999). Local regression and likelihood, volume 47. Springer New York.
#'
#' Franke, R. and Nielson, G. (1980). Smooth interpolation of large sets of scattered data.
#' International journal for numerical methods in engineering, 15(11):1691-1704.
#'
#' @seealso MGWRSAR, summary_mgwrsar, plot_mgwrsar, predict_mgwrsar
#' @examples
#' \donttest{
#' library(mgwrsar)
#' ## loading data example
#' data(mydata)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' mytab<-bandwidths_mgwrsar(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#' fixed_vars=c('Intercept','X1'),Models=c('GWR','MGWR'),Kernels=c('bisq','gauss'),
#' control=list(NN=300,adaptive=TRUE),control_search=list())
#'
#' names(mytab)
#' names(mytab[['GWR_bisq_adaptive']])
#'
#' mytab[['GWR_bisq_adaptive']]$config_model
#' mytab[['GWR_bisq_adaptive']]$CV
#' summary(mytab[['GWR_bisq_adaptive']]$model$Betav)
#'
#' mybestmodel=mytab[['GWR_gauss_adaptive']]$model
#' plot_mgwrsar(mybestmodel,type='B_coef',var='X2')
#' }
bandwidths_mgwrsar <- function(formula, data,coord, fixed_vars='Intercept',Models='GWR',Kernels_candidates='bisq',control=list(),control_search=list()){
  set.seed(123)
  if(sum(duplicated(coord))>0) {
    coord<-jitter(coord,0.01)
    warning('coords have been jittered because there is some duplicated location.',immediate. = TRUE)
  }
  ptmb<-proc.time()
  ### OLS
  if('OLS' %in% Models) Models<-Models[-which('OLS'==Models)]
  ### init
  n=nrow(data)
  k=length(attr(terms(as.formula(formula)),'variables'))-1
  con=list(searchB=T,Z=NULL,W=NULL,kernel_w=NULL,h_w=NULL,adaptive=T,Method='2SLS',TIME=FALSE,decay=0,Type='GD',isfgcv=TRUE,isgcv=FALSE,LocalInst='L5',Lambdacor=FALSE,NN=min(2000,n),Lambda=0, Lambdaj=rep(0,n),verbose=FALSE,remove_local_outlier=FALSE,outv=0.01,doMC=FALSE,ncore=1,Wh=NULL,SE=FALSE,TP=NULL,kWtp=16,KernelTP='sheppard',nstop=NULL,nneg=8,Wtp=NULL,tp_rmse=2)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
  for(i in 1:length(con))
  {
    assign(names(con)[i],con[[i]])
  }
  con_S=list(adaptive_W=F,kernel_w='rectangle',search_W=FALSE,Penalized=TRUE,n_searchW=1,verbose=TRUE)
  nmsC <- names(con_S)
  con_S[(namc <- names(control_search))] <- control_search
  if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in con_Strol: ", paste(noNms, collapse = ", "))
  for(i in 1:length(con_S))
  {
    assign(names(con_S)[i],con_S[[i]])
  }
  e_search=list()
  e_search$kernel_w=kernel_w
  e_search$search_W=search_W
  e_search$search_adaptive=adaptive_W
  e_search$n_searchW=n_searchW
  e_search$verbose=verbose
  e_search$tolerance<-tolerance<-0.0001
  e_search$Hp=NULL ## A voir

  if(Type!='GD') stop("bandwidths_mgwrsar is design only for 2D spatial kernel, i.e. Type='GD';\n
  For General Product Kernels, see bandwidths_mgwrsar_GPK() function")
  if(is.null(TP)) TP=1:n
  ######################

  DGPTAB<-init_DGPTAB()
  for(Model in Models) {
    if(verbose) cat('##### ', Model, ' #####\n')
    count=1
    stage1=prep_d(coord,NN,TP)
    con$indexG=stage1$indexG
    con$dists=stage1$dists
    for(kernels in Kernels_candidates) {
      if(verbose) cat('##### ', kernels, ' adaptive=',adaptive,' #####\n')
        if(!adaptive){
          lower=quantile(as.numeric(con$dists[,-1])[as.numeric(con$dists[,-1])>0],0.001)
          upper=max(as.numeric(con$dists))} else {
            if(kernels!='gauss') lower=2* ncol(attr(terms(as.formula(formula)),"factors")) else lower=3
            upper=NN-1
          }
        MyModel=paste0(Model,'_',kernels,ifelse(adaptive,'_adaptive',''))
        DGPTAB[[MyModel]]<-fb(formula,data,coord,fixed_vars,Model,con,kernels,e_search,lower,upper,tolerance)
        count=count+1
        if(verbose) cat('\n')
      }

  }
  DGPTAB[[MyModel]]$ctime=(proc.time()-ptmb)[3]
  DGPTAB[as.numeric(which(sapply(DGPTAB,length)>0))]
}
