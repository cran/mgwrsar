#' A bootstrap test for testing nullity of all Betas for mgwrsar class model,
#' @usage mgwrsar_bootstrap_test_all(model,B=100,ncore=1,type='standard')
#' @param model A mgwrsar model
#' @param B number of bootstrap replications, default 100
#' @param ncore number of cores for parallelization.
#' @param type type of boostrap ('spatial','wild','random')
#' @seealso  mgwrsar_bootstrap_test
#' @return a matrix with statistical test values and p ratios
mgwrsar_bootstrap_test_all <-
function(model,B=100,ncore=1,type='standard'){
x1=model
k=ncol(model@X)
names_x=colnames(model@X)
names_x<-names_x[-which(names_x=='Intercept')]
res<-matrix(nrow = length(names_x), ncol = 2)
for(i in 1:(k-1)){
name_x<-names_x[i]

x0<-MGWRSAR(formula=update(x1@formula,paste('~.-',name_x,sep='')),data=x1@data,coords=x1@coords,fixed_vars=x1@fixed_vars,kernels=x1@kernels,H=x1@H,Model=x1@Model,control=list(Method=x1@Method,W=x1@W,isgcv=FALSE,SE=TRUE))
cat(name_x,' ')
res[i,]<-unlist(mgwrsar_bootstrap_test(x0,x1,B=B,ncore=ncore,type=type,eps='H0',df='H0',focal='median',D=NULL))
}
colnames(res)<-c('Pratio','T')
as.data.frame(res)
}
