#' update_opt_known
#' to be documented
#' @usage update_opt_known(env = parent.frame())
#' @param model  to be documented
#' @noRd
#' @return to be documented
update_opt_known <-function(env = parent.frame()){with(env,{
  controlv=control
  data$e0k<-data$e0+BETA[,k]*X[,k]
  myformula_bk=as.formula(paste0('e0k~-1+',k))
  if(controlv$Type == "GDT"){
    if(opt[k]<max_dist & opt_t[k]<max_dist_t) model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(opt[k], opt_t[k]), Model = 'GWR',control=controlv)
    else if(opt[k]==max_dist & opt_t[k]==max_dist_t) model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=NULL, Model = 'OLS',control=controlv)
    else if(opt[k]<max_dist & opt_t[k]==max_dist_t) {
      if(!exists('controlvd', inherits = FALSE)) {
        controlvd=modifyList(controlv, list(dists=NULL,indexG=NULL,Type = 'GD',adaptive=controlv$adaptive[1]))
      }
      model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels,H=c(opt[k],NULL), Model = 'GWR',control=controlvd)
    }
    else if(opt[k]==max_dist & opt_t[k]<max_dist_t) {
      if(!exists('controlvt', inherits = FALSE)) {
        controlvt=modifyList(controlv, list(dists=NULL,indexG=NULL,Type = 'T',adaptive=F))
      }
      model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=NULL,kernels=kernels[2],H=opt_t[k], Model = 'GWR',control=controlvt)
    }
  } else {
    if((opt[k]>=n-2 & controlv$adaptive[1]) | (opt[k]>=max_dist & !controlv$adaptive[1]))  model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=k,kernels=kernels,H=c(opt[k],NULL), Model = 'OLS',control=controlv) else  model_k<-MGWRSAR(formula = myformula_bk, data = data,coords=coords, fixed_vars=k,kernels=kernels,H=c(opt[k],NULL), Model = 'GWR',control=controlv)
  }

  e0=residuals(model_k)
  betav<-model_k@Betav
  isol=is.na(betav)
  e0[isol] <- data$e0[isol]
  if(get_AIC ) {
    TSik[!isol,k]<-model_k@TS
    Sk<-model_k@Shat
  }
})
}




