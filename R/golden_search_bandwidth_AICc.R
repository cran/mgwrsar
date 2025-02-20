#' golden_search_bandwidth
#' to be documented
#' @usage golden_search_bandwidth(formula,H2=NULL,data, coords, fixed_vars,
#' kernels, Model, control,lower.bound, upper.bound,tolerance=0.000001)
#' @param formula to be documented
#' @param H2 to be documented
#' @param data to be documented
#' @param coords to be documented
#' @param fixed_vars to be documented
#' @param kernels to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param lower.bound to be documented
#' @param upper.bound to be documented
#' @param tolerance to be documented
#' @return a list(minimum=res,objective=objective,model=model).
golden_search_bandwidth<-function(formula,H2=NULL,data, coords, fixed_vars, kernels, Model, control,lower.bound, upper.bound,tolerance=0.000001)
{
  ptm=proc.time()
  set.seed(123)
  while(sum(duplicated(coords))>0) {
    coords<-jitter(coords,0.001)
    #warning('coords have been jittered because there is some duplicated location.')
  }

  mycall_golden_search_bandwidth <- match.call()
  if(control$verbose) cat('.')
  adaptive=control$adaptive[1]
  golden.ratio = 2/(sqrt(5) + 1)
  ### Use the golden ratio to set the initial test points
  x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  if(adaptive) {
    x1=floor(x1)
    x2=ceiling(x2)
  }
  ## distance computation
  if(is.null(control$dists)){
    if(is.null(control$NN)) control$NN=nrow(data)
    if(is.null(control$TP)) control$TP=1:nrow(data)
    if (length(kernels) > 1) S = as.matrix(cbind(coords, control$Z)) else S = as.matrix(coords)
    stage1=prep_d(S,control$NN,control$TP)
    control$indexG=stage1$indexG
    control$dists=stage1$dists
  }
  ### On evalube cv_h
  H<-c(x1,H2)
  f1 = AICc_CV(H,formula, data ,coords, fixed_vars,kernels, Model,control)

  H<-c(x2,H2)
  f2 = AICc_CV(H,formula, data ,coords, fixed_vars,kernels, Model,control)

  iteration = 0
  if(adaptive) tolerance=1
  x1_p=1
  x2_p=2
  while (abs(upper.bound - lower.bound) > tolerance & ((adaptive & abs(x1-x2)>1 & (x1_p!=x1 | x2_p!=x2) ) | !adaptive))
  {
    if(control$verbose) cat('.')
    if(control$verbose)  cat('f1=',f1[1],' f2=',f2[1],' x1=',x1[1],' x2=',x2[1],'\n')
    x1_p=x1
    x2_p=x2
    #if(iteration==4) browser()
    iteration = iteration + 1
    ## f2>f1
    if(f2 > f1 & abs(upper.bound - lower.bound) > tolerance){
      upper.bound = x2
      x2 = x1
      f2=f1
      x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
      if(adaptive) x1=floor(x1)
      H<-c(x1,H2)
      f1 = AICc_CV(H,formula, data ,coords, fixed_vars,kernels, Model,control)
    }
    ## f2<=f1
    if(f2 <= f1 & abs(upper.bound - lower.bound) > tolerance){
      #if(control$verbose)  cat('f1=',f1[1],' f2=',f2[1],' x1=',x1[1],' x2=',x2[1],'\n')
      lower.bound = x1
      x1 = x2
      f1=f2
      x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
      if(adaptive) x2=ceiling(x2)
      H<-c(x2,H2)
      f2 = AICc_CV(H,formula, data ,coords, fixed_vars,kernels, Model,control)
    }
  }
  res=(lower.bound+upper.bound)/2
  if(adaptive ) {
    if( f1<f2) {
      res=x1
      objective=f1
    } else {
      res=x2
      objective=f2
    }
  } else  objective=AICc_CV(H,formula, data ,coords, fixed_vars,kernels, Model,control)
  ctime=(proc.time()-ptm)[3]
  model<-MGWRSAR(formula, data,coords, fixed_vars,kernels,H=c(res,H2), Model = Model,control=control)
  model@mycall[['formula']]<-mycall_golden_search_bandwidth[['formula']]
  model@mycall[['data']]<-mycall_golden_search_bandwidth[['data']]
  model@mycall[['coords']]<-mycall_golden_search_bandwidth[['coords']]
  model@mycall[['kernels']]<-mycall_golden_search_bandwidth[['kernels']]
  model@mycall[['control']]<-mycall_golden_search_bandwidth[['control']]
  list(minimum=res,objective=objective,model=model,ctime=ctime)
}
