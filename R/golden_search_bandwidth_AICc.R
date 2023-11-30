#' golden_search_bandwidth
#' to be documented
#' @usage golden_search_bandwidth_AICc(formula, data, coords, fixed_vars,
#' kernels, Model, control,lower.bound, upper.bound,tolerance=1)
#' @param formula to be documented
#' @param data to be documented
#' @param coords to be documented
#' @param fixed_vars to be documented
#' @param kernels to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param lower.bound to be documented
#' @param upper.bound to be documented
#' @param tolerance to be documented
#' @noRd
#' @return a list(minimum=res,objective=objective,model=model).
golden_search_bandwidth_AICc<-function(formula, data, coords, fixed_vars, kernels, Model, control,lower.bound, upper.bound,tolerance=1)
{
  cat('\n .')
  adaptive=control$adaptive
  golden.ratio = 2/(sqrt(5) + 1)
  ### Use the golden ratio to set the initial test points
  x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  if(adaptive) {
    x1=floor(x1)
    x2=ceiling(x2)
  }
  ### On evalube cv_h
  H<-x1
  f1 = AICc(H,formula, data ,coords, fixed_vars,kernels, Model,control)

  H<-x2
  f2 = AICc(H,formula, data ,coords, fixed_vars,kernels, Model,control)

  iteration = 0
  if(adaptive) tolerance=1
  x1_p=1
  x2_p=2
  while (abs(upper.bound - lower.bound) > tolerance & ((adaptive & abs(x1-x2)>1 & (x1_p!=x1 | x2_p!=x2) ) | !adaptive))
  {
    cat('.')
    #cat('f1=',f1[1],' f2=',f2[1],' x1=',x1[1],' x2=',x2[1],'\n')
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
      H<-x1
      f1 = AICc(H,formula, data ,coords, fixed_vars,kernels, Model,control)
    }
    ## f2<=f1
    if(f2 <= f1 & abs(upper.bound - lower.bound) > tolerance){
      lower.bound = x1
      x1 = x2
      f1=f2
      x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
      if(adaptive) x2=ceiling(x2)
      H<-x2
      f2 = AICc(H,formula, data ,coords, fixed_vars,kernels, Model,control)
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
  } else  objective=AICc(H,formula, data ,coords, fixed_vars,kernels, Model,control)
  model<-MGWRSAR(formula, data,coords, fixed_vars,kernels,H=res, Model = Model,control=control)
  list(minimum=res,objective=objective,model=model)
}
