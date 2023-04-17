#' fb
#' to be documented
#' @usage fb(formula,data,coord,fixed_vars,Model,control,
#' kernels,e_search,lower,upper,tolerance)
#' @param formula to be documented
#' @param data to be documented
#' @param coord to be documented
#' @param fixed_vars to be documented
#' @param Model to be documented
#' @param control to be documented
#' @param kernel to be documented
#' @param e_search to be documented
#' @param lower to be documented
#' @param upper to be documented
#' @param tolerance  to be documented
#' @noRd
#' @return to be documented
fb <-
  function(formula,data,coord,fixed_vars,Model,control,kernels,e_search,lower,upper,tolerance){
    control$isgcv<-TRUE
    starting=10^100

    # if(e_search$search_W){
    # voir old fb
    # }
    if (e_search$search_W  & !(Model %in% c('GWR','OLS','MGWR'))) {
     if(e_search$verbose) cat("\n########\nSearch W stage\n########")
      if(e_search$search_adaptive) {
        upperW=50
        lowerW=1
        } else {
          upperW=max(control$dists[,50])
          lowerW=median(control$dists[,2])
        }
      #cat('\n########\nSearch h for GWR model \n########\n')
      if(Model!='SAR'){
      opt=golden_search_bandwidth(Hp=NULL,kernel_w=NULL,search_adaptive=e_search$search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=ifelse(is.null(fixed_vars),'GWR','MGWR'), control=control,lower.bound=lower, upper.bound=upper,tolerance=tolerance)
      h=opt$minimum
      #cat(' objective = ',opt$objective,' minimum = ',opt$minimum,'\n')
      } else h=NULL
      optW=golden_search_bandwidth(Hp=10*h,kernel_w=e_search$kernel_w,search_adaptive=e_search$search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control,lower.bound=lowerW, upper.bound=upperW,tolerance=tolerance)
      if(e_search$search_adaptive) NNN=optW$minimum+2 else NNN=500
      control$W=kernel_matW(H=optW$minimum,kernels=e_search$kernel_w,coord_i=coord,NN=NNN,adaptive=e_search$search_adaptive,diagnull=TRUE)
      if(e_search$verbose) cat('\n W : kernel =',e_search$kernel_w,ifelse(e_search$search_adaptive,'adaptive',''),': objective = ',optW$objective,' minimum = ',optW$minimum,'\n\n')
    }

    if(Model!='SAR'){
      if(e_search$verbose) cat('\n########\nSearch bandwidth stage\n########\n')
      control$isgcv=TRUE
      #stage1
      if(control$adaptive) uppermin=0.1*control$NN else  uppermin=upper
 opt=golden_search_bandwidth(Hp,kernel_w=NULL,search_adaptive=e_search$search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control,lower.bound=lower, upper.bound=uppermin,tolerance=tolerance)
 if(opt$minimum>uppermin*0.95 & control$adaptive){
   if(e_search$verbose) cat('\n ...try larger suppport\n')
   opt=golden_search_bandwidth(Hp,kernel_w=NULL,search_adaptive=e_search$search_adaptive,formula=formula, data=data, coord=coord, fixed_vars=fixed_vars, kernels=kernels, Model=Model, control=control,lower.bound=lower, upper.bound=upper,tolerance=tolerance)
   if(opt$minimum>uppermin*0.95 & control$adaptive) cat('\n Border solution !!! Try to increase NN')
 }
 if(e_search$verbose) cat('kernel =',kernels,ifelse(control$adaptive,'adaptive',''),' objective = ',opt$objective,' minimum = ',opt$minimum,'\n')
      config_model<-list(Model=Model,kernels=kernels[1],adaptive=control$adaptive,H=opt$minimum,kernel_w=e_search$kernel_w,h_w=e_search$h_w)
    } else config_model<-list(Model=Model,kernels=NULL,adaptive=NULL,H=NULL,kernel_w=e_search$kernel_w,h_w=e_search$h_w)
    fill_DGPTAB(formula,data,coord,fixed_vars,Model,control,opt,e_search$search_W,config_model)
  }
