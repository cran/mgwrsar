#' init_param_tds
#' to be documented
#' @usage init_param_tds(env = parent.frame())
#' @param env an environment
#' @return an environment
#' @noRd
init_param_tds<-function(env = parent.frame()){
  with(env,{
  if(!(Model %in% c('tds_mgwr','atds_mgwr','atds_gwr'))) stop('Only atds_gwr, tds_mgwr and  atds_mgwr Model can be estimated using Top Down Scale approach in this release.')
  n_time<-m<-n<-nrow(data)
  # if(!is.null(control$TP)) {
  #   if(length(control$TP)!=n)
  #     cat('tds_mgwr is not designed for Target Point, automatic switch to control$TP=1:n ')
  #   #control$TP=1:n
  # } else control$TP=1:n
  ### param tds
  type='proportional'
  if(is.null(control_tds$init_model)) control_tds$init_model='OLS'
  if(is.null(control_tds$nns)) control_tds$nns=30
  if(is.null(control_tds$get_AIC) ) control_tds$get_AIC=FALSE
  if(Model=='atds_mgwr')  control_tds$get_AIC=TRUE
  if(is.null(control_tds$doMC)) control_tds$doMC=FALSE
  if(is.null(control_tds$ncore)) control_tds$ncore=1
  if(is.null(control_tds$tol)) control_tds$tol=0.001
  if(is.null(control_tds$nrounds)) control_tds$nrounds=3
  if(is.null(control_tds$verbose)) control_tds$verbose=FALSE
  if(is.null(control_tds$V)) V<-control_tds$V<-NULL
  if(is.null(control_tds$min_dist)) min_dist=NULL
  if(is.null(control_tds$first_nn)) control_tds$first_nn=n
  if(is.null(control_tds$BETA)) control_tds$BETA=NULL
  if(is.null(control_tds$TRUEBETA)) control_tds$TRUEBETA=NULL
  if(is.null(control_tds$H)) H<-control_tds$H<-NULL
  if(is.null(control_tds$browser)) control_tds$browser=0
  reassign_control(control_tds)
  if(!('model_stage1' %in% ls())) model_stage1=NULL
  control$get_ts=TRUE
  control$get_s=get_AIC
  if(control$doMC & control_tds$doMC) {
    cat('WARNING: nested parallel loops not supported, control_tds$doMC switch off.')
    control_tds$doMC=FALSE
  }
  if(is.null(control$isgcv)) control$isgcv=FALSE
  if(is.null(control$Type)) control$Type='GD'
  if(is.null(control$NN)) control$NN=min(3002,n)
  if(is.null(control$TP)) control$TP=1:n
  if(is.null(control$alpha)) control$alpha=1
  if(is.null(control$theta)) control$theta=1
  if(is.null(control$doMC)) control$doMC=FALSE
  if(control_tds$doMC & control$doMC & control_tds$ncore>1) control$doMC<- FALSE


  while(sum(duplicated(coords))>0) {
    coords<-jitter(coords,0.001)
    #warning('coords have been jittered because there is some duplicated location.')
  }

  if(!is.null(control$Z)) coords_in=as.matrix(cbind(coords,control$Z)) else coords_in=coords
  if(!('indexG' %in% names(control)) & is.null(model_stage1)) { G<-prep_d(coords_in,control$NN,control$TP)
  control$indexG=G$indexG
  control$dists=G$dists
  } else if(!is.null(model_stage1)){
    if(is.null(model_stage1@G)) G<-prep_d(coords_in,control$NN,control$TP) else G<-model_stage1@G
    control$indexG=G$indexG
    control$dists=G$dists
  }
  #browser()
  if(control$adaptive[1]) max_dist=n else max_dist=max(control$dists[[1]][,ncol(control$dists[[1]])])


  ### model
  mycall <- match.call()
  mf <- model.frame(formula, data)
  mt <- attr(x = mf, which = "terms")
  Y <- model.extract(mf, "response")
  X = model.matrix(object = mt, data = mf)
  colnames(X)<-clean_colnames(X)
  if(!is.null(fixed_vars)) fixed_vars=intersect(fixed_vars,colnames(X))
  if(ncol(X)>1) "tds_gwr is only designed for univariate regression. With multivariate GWR, you cannot be sure of achieving a global optimum. The results may show better in-sample RMSE or AICc compared to GWR, tds_mgwr, or atds_mgwr, but it must be compared to these models using cross-validation."


  formula<-formula(lm(formula,data))
  if(!('model_stage1' %in% ls())) model_stage1=NULL
  if(!is.null(model_stage1)) {
    S<-model_stage1@Shat
    if(!(model_stage1@Model %in% c('OLS','GWR'))) V<-control_tds$V<-model_stage1@V
    BETA=model_stage1@Betav
    AICc=model_stage1@AICc
  }
  I<-Diagonal(n)
  namesX =colnames(X)
  if(colnames(X)[1]== "(Intercept)") colnames(X)[1]<-namesX[1]<-'Intercept'
  data$Intercept=1
  K=length(namesX)
  if(kernels[1]=='gauss') control_tds$minv<-minv<-2 else if(is.null(control_tds$minv)) control_tds$minv<-minv<-K+1
  if(!('BETA' %in% ls())) BETA=NULL
  if(!('TRUEBETA' %in% ls())) TRUEBETA=NULL
  varying<-setdiff(namesX,fixed_vars)
  if(is.null(BETA)){
    if(Model %in% c('tds_mgwr','atds_mgwr')) starting_model_tds_mgwr() else if(Model=='atds_gwr')  starting_model_atds_gwr()
  }
  myformula_b=as.formula(paste0('e0~-1+',paste0(varying,collapse = '+')))
  if(length(fixed_vars)>0) formula_constant=as.formula(paste0('e0~-1+',paste0(fixed_vars,collapse = '+')))

  #if(length(H)>0) names(H)<-namesX

  # algo History
  HBETA=list()
  HBETA[[1]]<-BETA
  if(get_AIC){
  HAICc<-c()
    HAICc<-c(HAICc,AICc)
    HTS=list()
  }

  HOPT=rep(NA,K)
  names(HOPT)<-namesX
  })
}

