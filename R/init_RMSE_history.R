# RMSE History for monte carlo
#' init_RMSE_history
#' to be documented
#' @usage built_Vseq(d,h)
#' @param env an environment
#' @return an environment.
#' @noRd
#'
init_RMSE_history<-function(env=parent.frame()){
  with(env,{
    if(!is.null(TRUEBETA)){
    HRMSE=matrix(NA,nrow=length(V)*2+1,ncol=K+5)
    for(k in 1:K) HRMSE[1,k]=sqrt(mean((TRUEBETA[,k]-BETA[,k])^2))
    HRMSE[1,K+1]<-mean(HRMSE[1,1:K])
    HRMSE[1,K+2]<-m
    if(!is.null(model_lm0)) HRMSE[1,K+3]<-model_lm0@AIC
    HRMSE[1,K+4]<-sqrt(mean(data$e0^2))
  }
  })
}
