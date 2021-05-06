#' Select optimal kernel and bandwidth from a list of models, kernels and bandwidth candidates.
#'
#' Given a lm formula and a dataframe with coordinates, function bandwidths_mgwrsar optimizes the choice of
#' a bandwidth value for each of the chosen models and kernel types using a leave-one-out cross validation criteria.
#' A cross validated criteria is also used for selecting the best kernel type for a given model.
#'
#' @usage bandwidths_mgwrsar(formula, data,coord,
#' fixed_vars='Intercept',Models='GWR',Kernels='bisq',control=list(),control_search=list())
#'
#' @param formula  a formula.
#' @param data a dataframe or a spatial dataframe (sp package).
#' @param coord a dataframe or a matrix with coordinates, not required if data is a spatial dataframe, default NULL.
#' @param fixed_vars a vector with the names of spatially constant coefficient. For mixed model, if NULL, the default
#' #' is set to 'Intercept'.
#' @param Models character containing the type of model: Possible values are "OLS",
#' "SAR", "GWR" (default), "MGWR" , "MGWRSAR_0_0_kv","MGWRSAR_1_0_kv",
#' "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0".
#' @param Kernels a vector with the names of kernel type.
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
#' @seealso MGWRSAR, summary_mgwrsar, plot_mgwrsar, predict_mgwrsar, kernelW_C
#'
#' @examples
#' \donttest{
#' library(mgwrsar)
#' data(mydata)
#' coord=as.matrix(mydata[,c("x_lat","y_lon")])
#' W=KNN(coord,8)
#' ######################
#' #### Finding bandwith by hand
#' #####################
#'
#' ### kernel only space
#' model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#' fixed_vars=NULL,kernels=c('bisq_knn'),H=50,
#' Model = 'GWR', control=list(isgcv=FALSE,minv=1))
#' cat('SSR =')
#' summary_mgwrsar(model_GWR)
#'
#' myCV<-function(H){model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,
#' coord=coord, fixed_vars=NULL,kernels=c('gauss_adapt'),H,
#'  Model = 'GWR',control=list(isgcv=TRUE))
#' model_GWR$SSR
#' }
#'
#' res=optimize(myCV,upper=500,lower=10)
#' res
#'
#' ## model with optimal bandwith with adaptative gaussian kernel
#'
#' model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#' fixed_vars=NULL,kernels=c ('gauss_adapt'),H=ceiling(res$minimum),
#' Model = 'GWR',control=list(isgcv=FALSE))
#' summary_mgwrsar(model_GWR)
#'
#' ######################
#' #### finding the bandwidths using bandwidths_mgwrsar
#' #####################
#'
#'
#' mytab<-bandwidths_mgwrsar(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#' fixed_vars='Intercept',Models=c('GWR','MGWR'),Kernels=c('bisq_knn','gauss_adapt','gauss'),
#' control=list(),control_search=list (lower_d=8,lower_c=0.03,upper_c=0.65))
#'
#' names(mytab)
#' names(mytab[['GWR']])
#'
#' mytab[['GWR']]$config_model
#' mytab[['GWR']]$CV
#' summary(mytab[['GWR']]$model$Betav)
#'
#'
#' mytab[['GWR_2']]$config_model
#' mytab[['GWR_2']]$CV
#' summary(mytab[['GWR_2']]$model$Betav)
#'
#' mytab[['GWR_3']]$config_model
#' mytab[['GWR_3']]$CV
#' summary(mytab[['GWR_3']]$model$Betav)
#'
#' mytab[['MGWR']]$config_model
#' mytab[['MGWR']]$CV
#' mytab[['MGWR']]$model$Betac
#' summary(mytab[['MGWR']]$model$Betav)
#'
#' mytab[['MGWR_2']]$config_model
#' mytab[['MGWR_2']]$CV
#' mytab[['MGWR_2']]$model$Betac
#' summary(mytab[['MGWR_2']]$model$Betav)
#'
#' mytab[['MGWR_3']]$config_model
#' mytab[['MGWR_3']]$CV
#' mytab[['MGWR_3']]$model$Betac
#' summary(mytab[['MGWR_3']]$model$Betav)
#'
#' mytab2<-bandwidths_mgwrsar(formula = 'Y_gwr~X1+X2+X3', data = mydata,coord=coord,
#'  fixed_vars='Intercept',Models=c('MGWRSAR_0_kc_kv'),Kernels=c('gauss_adapt'),
#'  control=list(),control_search=list(search_W=TRUE,kernels_w=c('bisq','gauss_adapt')))
#'
#' mytab2[['MGWRSAR_0_kc_kv']]$config_model
#' mytab2[['MGWRSAR_0_kc_kv']]$CV
#' mytab2[['MGWRSAR_0_kc_kv']]$model$Betac
#' summary(mytab2[['MGWRSAR_0_kc_kv']]$model$Betav)
#'
#' }
bandwidths_mgwrsar <-
function(formula, data,coord, fixed_vars='Intercept',Models='GWR',Kernels='bisq',control=list(),control_search=list()){
### OLS
if('OLS' %in% Models) Models<-Models[-which('OLS'==Models)]
### init
n=nrow(data)
k=length(attr(terms(as.formula(formula)),'variables'))-1

con=list(Z=NULL,W=NULL,kernel_w=NULL,h_w=NULL,Method='2SLS',TIME=FALSE,decay=0,minv=30,Type='GD',isfgcv=TRUE,isgcv=FALSE,LocalInst='L5',Lambdacor=FALSE,maxknn=500,NmaxDist=5000,Lambda=0, Lambdaj=rep(0,n),verbose=FALSE,remove_local_outlier=FALSE,outv=0.01,doMC=FALSE,ncore=1,Wh=NULL,SE=FALSE,xratiomin=10e-10)

nmsC <- names(con)
con[(namc <- names(control))] <- control
if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in control: ", paste(noNms, collapse = ", "))
for(i in 1:length(con))
{
assign(names(con)[i],con[[i]],envir =parent.frame())
}


ds<-sample(1:n,30)
D1=as.matrix(dist(coord))
P01D=quantile(D1,prob=0.01)
P99D=quantile(D1,prob=0.999)
P005D=quantile(D1,prob=0.005)
P20D=quantile(D1,prob=0.20)
con_S=list(kernels_w=c('gauss_adapt','bisq_knn'),search_W=FALSE,lower_c=P01D,upper_c=P99D,lower_d=2*k+1,lower_cW=P005D,lower_dW=2,upper_dW=2,Penalized=TRUE,n_searchW=1)
nmsC <- names(con_S)
con_S[(namc <- names(control_search))] <- control_search
if (length(noNms <- namc[!namc %in% nmsC]))  warning("unknown names in con_Strol: ", paste(noNms, collapse = ", "))
for(i in 1:length(con_S))
{
assign(names(con_S)[i],con_S[[i]],envir =parent.frame())
}
e_search=list()
e_search$kernels_w=kernels_w
e_search$lower_cw=lower_c
e_search$upper_cw=upper_c
e_search$lower_dw=lower_d
e_search$lower_dW=lower_dW
e_search$lower_cW=lower_cW
e_search$upper_dW=upper_dW
e_search$Penalized=Penalized
e_search$search_W=search_W
e_search$n=n
e_search$kernel_w=('unknow')
e_search$h_w=0
e_search$n_searchW=n_searchW

######################

DGPTAB<-init_DGPTAB()
for(Model in Models) {
cat('##### ', Model, ' #####\n')
	count=1
	for(kernel in Kernels) {
		cat('##### ', kernel, ' #####\n')
		if(count==1) MyModel=Model else MyModel=paste(Model,count,sep='_')
		DGPTAB[[MyModel]]<-fb(formula,data,coord,fixed_vars,Model,con,kernel,e_search)
		count=count+1
		cat('\n')
	}
}
DGPTAB[as.numeric(which(sapply(DGPTAB,length)>0))]
}
