#' A bootstrap test for Betas for mgwrsar class model.
#' @usage mgwrsar_bootstrap_test(x0,x1,B=100,domc=FALSE,ncore=1,
#' type='standard',eps='H1',df='H1',focal='median',D=NULL)
#' @param x0  The H0 mgwrsar model
#' @param x1  The H1 mgwrsar model
#' @param B  number of bootstrap repetitions, default 100
#' @param domc  If TRUE, doParallel parallelization
#' @param ncore  number of cores
#' @param type  type of bootstap : 'wild','Rademacher','spatial' or 'standard' (default)
#' @param eps  Hypothesis under wich residuals are simulated,  'H0' or 'H1' (default)
#' @param df  Hypothesis under wich degree of freedom is estimated.
#' @param focal  see sample_stat help
#' @param D  A matrix of distance
#' @seealso  mgwrsar_bootstrap_test_all
#' @return The value of the statictics test and a p ratio.
#' @examples
#' \donttest{
#'  data(data_mgwrsar)
#'  coord=as.matrix(mydata[,c("x_lat","y_lon")])
#'  model_GWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata, coord=coord,
#'   fixed_vars=NULL,kernels=c('gauss_adapt'), H=20,
#'   Model = 'GWR',control=list(SE=TRUE))
#'
#'  model_MGWR<-MGWRSAR(formula = 'Y_gwr~X1+X2+X3', data = mydata, coord=coord,
#'  fixed_vars='Intercept',kernels=c('gauss_adapt'), H=20,
#'  Model = 'MGWR',control=list(SE=TRUE))
#'
#'  mgwrsar_bootstrap_test(model_MGWR,model_GWR,domc=FALSE,B=30)
#' }
mgwrsar_bootstrap_test <-
function(x0,x1,B=100,domc=FALSE,ncore=1,type='standard',eps='H1',df='H1',focal='median',D=NULL) {
if(is.null(D) & type=='spatial') stop("When type='spatial' ou have to provide D a nxn distance matrix")
	ptm<-proc.time()
	n=length(x0$Y)
    T<-Tf(x0,x1)
    if(df=='H1' & eps=='H1') eps1<-(x1$residuals/sqrt(n/x1$edf)) -mean(x1$residuals/sqrt(n/x1$edf)) else if(df=='H1'& eps=='H0') eps1<-(x0$residuals/sqrt(n/x1$edf)) -mean(x0$residuals/sqrt(n/x1$edf))  else if(df=='H0'& eps=='H0') eps1<-(x0$residuals/sqrt(n/x0$edf)) -mean(x0$residuals/sqrt(n/x0$edf))  else eps1<-(x1$residuals/sqrt(n/x0$edf)) -mean(x1$residuals/sqrt(n/x0$edf))

    if(domc==FALSE){
    T_star<-numeric(B)
    for(i in 1:B){
        if(type=='wild') {
        eps2=rwild(eps1,"Rademacher")
    } else if(type=='spatial') {ind<-sample_spat(D,x0$H,focal)
    eps2<-eps1[ind]} else {eps2<-as.matrix(sample(eps1, replace=TRUE))}
        Y_star=x0$fit+eps2
        simWS=x0$data
   		  simWS[, all.vars(as.formula(x0$formula))[1]]=Y_star

   x_star0<-MGWRSAR(formula=x0$formula,data=simWS,coord=x0$coord,fixed_vars=x0$fixed_vars,kernels=x0$kernel,H=x0$H,Model=x0$Model,control=list(Method=x0$Method,W=x0$W,isgcv=x0$isgcv,SE=TRUE))
   x_star1<-MGWRSAR(formula=x1$formula,data=simWS,coord=x1$coord,fixed_vars=x1$fixed_vars,kernels=x1$kernel,H=x1$H,Model=x1$Model,control=list(Method=x1$Method,W=x1$W,isgcv=x1$isgcv,SE=TRUE))
        T_star[i]=Tf(x_star0,x_star1)
        cat('T_star[',i,'] ',T_star[i],' ')
    	}
    } else
    {
    registerDoParallel()
    T_star<-foreach(i=1:B,.combine="rbind",.inorder=FALSE)  %dopar% {
    	          if(type=='wild') {
                eps=rwild(eps1,"Rademacher")
    } else if(type=='spatial') {ind<-sample_spat(D,x0$H,focal)
    eps<-eps1[ind]} else {eps<-as.matrix(sample(eps1, replace=TRUE))}
        Y_star=x0$fit+eps
        simWS=x0$data
        simWS[, all.vars(as.formula(x0$formula))[1]]=Y_star

  		x_star0<-MGWRSAR(formula=x0$formula,data=simWS,coord=x0$coord,fixed_vars=x0$fixed_vars,kernels=x0$kernel,H=x0$H,Model=x0$Model,control=list(Method=x0$Method,W=x0$W,isgcv=x0$isgcv,SE=TRUE))
   		x_star1<-MGWRSAR(formula=x1$formula,data=simWS,coord=x1$coord,fixed_vars=x1$fixed_vars,kernels=x1$kernel,H=x1$H,Model=x1$Model,control=list(Method=x1$Method,W=x1$W,isgcv=x1$isgcv,SE=TRUE))
      Tf(x_star0,x_star1)
    	}

    }
    p_star=sum(T_star>T)/B
    diff<-proc.time()-ptm
    list(p_star=p_star,T=T)
}
