#' starting_model_atds_gwr
#' to be documented
#' @usage starting_model_atds_gwr(env = parent.frame())
#' @param env an environment
#' @return an environment.
#' @noRd
starting_model_atds_gwr<-function(env = parent.frame()){
with(env,{
#browser()
model_lm0<-MGWRSAR(formula = formula, data = data,coords=coords,fixed_vars=NULL, Model = 'OLS',H=NULL,kernels=NULL,control=list(get_s=TRUE))
BETA=matrix(coef(model_lm0)$Betac,byrow=TRUE,nrow=,nrow(data),ncol=K)
colnames(BETA)<-namesX
data$e0cv<-data$e0<-residuals(model_lm0)
idx_init<-idx<-1:n
drop_init=c()
n_time<-n
global_ts<-sum(model_lm0@TS)
S<-model_lm0@Shat
AICc<-n*log(sum(data$e0^2)/n)+n*log(2*pi)+n*(n+global_ts)/(n-2-global_ts)
})
}
