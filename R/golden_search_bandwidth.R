#' golden_search_bandwidth
#' to be documented
#' @usage golden_search_bandwidth(Hp,kernel_w,search_adaptive,formula,data,
#' coord,fixed_vars,kernels,Model,control,lower.bound, upper.bound,tolerance)
#' @param Hp to be documented
#' @param kernel_w to be documented
#' @param search_adaptive to be documented
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param kernels to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param lower.bound to be documented
#' @param upper.bound to be documented
#' @param tolerance to be documented
#' @noRd
#' @return a vector of weights.
golden_search_bandwidth<-function(Hp,kernel_w,search_adaptive,formula, data, coord, fixed_vars, kernels, Model, control,lower.bound, upper.bound,tolerance)
{
  cat('\n .')
  adaptive=control$adaptive
  golden.ratio = 2/(sqrt(5) + 1)
  if(!is.null(kernel_w)) adaptive=control$adaptive
  ### Use the golden ratio to set the initial test points
  x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  if(adaptive) {
    x1=floor(x1)
    x2=ceiling(x2)
  }
  ### On evalube cv_h
  if(is.null(kernel_w)){Hp<-H<-x1 } else {H<-x1 }

  f1 = cv_h(H,Hp,kernel_w=kernel_w,search_adaptive=search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control)

  if(is.null(kernel_w)){Hp<-H<-x2 } else {H<-x2 }
  f2 = cv_h(H,Hp,kernel_w=kernel_w,search_adaptive=search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control)

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
      if(is.null(kernel_w)){Hp<-H<-x1 } else {H<-x1 }
      f1 = cv_h(H,Hp,kernel_w=kernel_w,search_adaptive=search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control)
    }
    ## f2<=f1
    if(f2 <= f1 & abs(upper.bound - lower.bound) > tolerance){
      lower.bound = x1
      x1 = x2
      f1=f2
      x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
      if(adaptive) x2=ceiling(x2)
      if(is.null(kernel_w)){Hp<-H<-x2 } else {H<-x2 }
      f2 = cv_h(H,Hp,kernel_w=kernel_w,search_adaptive=search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control)
    }
  }
  res=(lower.bound+upper.bound)/2
  if(is.null(kernel_w)) Hp<-H
  if(adaptive ) {
    if( f1<f2) {
      res=x1
      objective=f1
    } else {
      res=x2
      objective=f2
    }
} else  objective=cv_h(H=res,Hp=Hp,kernel_w=kernel_w,search_adaptive=search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control)
  list(minimum=res,objective=objective)
}


