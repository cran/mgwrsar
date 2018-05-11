#' GWR
#' to be documented
#' @usage GWR(Y,XV,XX,YY,S,H, kernels, type = "GD", minv = 1, maxknn = 500,
#'  NmaxDist = 6000,SE=FALSE, isgcv, TIME, decay,interceptv=TRUE,W=NULL,
#'  betacor=FALSE,remove_local_outlier=FALSE,outv=0,doMC=FALSE,ncore=1,
#'  Wh=NULL,xratiomin=10e-10)
#' @param Y  A vector of response
#' @param XV  A matrix with covariates with stationnary parameters
#' @param XX   A matrix with all covariates (XC,XV)
#' @param YY  A vector of first stage Y for some models (see MGWRSAR code)
#' @param S  A matrix with variables used in kernel
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param type  Type of Genelarized kernel product ('GD' only spatial,'GDC' spatial + a categorical variable,'GDX' spatial + a continuous variable,'GDT' spatial + a time index, and other combinations 'GDXXC','GDTX',...)
#' @param minv  Minimum number of non null weight
#' @param NmaxDist  Maximum number of observation for computing dense weight matrix
#' @param maxknn  If n >NmaxDist how many column with dense weight matrix (max number of neighbours)
#' @param SE  If standard error are computed, default FALSE
#' @param isgcv  leave one out cross validation, default FALSE.
#' @param TIME  Use rigth truncated kernel for time index kernel
#' @param decay  time decay
#' @param interceptv  Intercept spatially varying, default FALSE
#' @param W  A weight matrix for spatial autocorrelation
#' @param betacor Do a tuncation of spatial autocorelation if absolute value larger than 1.
#' @param remove_local_outlier Remove local outlier
#' @param outv  A treshold for removing local outlier
#' @param doMC  doParallel parallelization
#' @param ncore  Number of cores for parallelization
#' @param Wh  A matrix of weights for local estimation
#' @param xratiomin  A treshold parameters for removing obs with not enough positive weigths for local regression
#' @return a list of object for MGWRSAR wrapper
#' @keywords internal
GWR <-
function(Y,XV,XX,YY,S,H, kernels, type = "GD", minv = 1, maxknn = 500, NmaxDist = 6000,SE=FALSE, isgcv, TIME, decay,interceptv=TRUE,W=NULL,betacor=FALSE,remove_local_outlier=FALSE,outv=0,doMC=FALSE,ncore=1,Wh=NULL,xratiomin=10e-10)
{
SEV=NULL
X=XV
n<-nrow(X)
if(SE) {tS<-0;edf=n;}
m<-ncol(X)
X=as.matrix(X)
if(!is.null(W)) PhWy=PhWY_C(as.matrix(YY),as.matrix(XX),W,rep(1,n))

if(doMC==FALSE){
Betav <- matrix(0,nrow=n, ncol=ifelse(is.null(W),m,m+1))
if(SE)  SEV <- matrix(0,nrow=n, ncol=ifelse(is.null(W),m,m+1))
for (i in 1:n){
if(is.null(Wh)) {
w.i<-GPKj(i-1,Y,X,S, H, kernels, type, minv = minv, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay) } else w.i<-Wh[i,]
index=which(w.i>0) #index=which((w.i[w.i>0]*length(w.i))>wimin)
Wd<-sqrt(w.i[index])
YYw<-Wd*Y[index]
keep=1:ncol(X)
if(m>1) {
toremove=unique(c(which(apply(X[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))),which((apply(Wd*X[index,],2,function(x) sum(abs(x)))/apply(XV[index,],2,function(x) sum(abs(x))))<xratiomin)))
if(interceptv) toremove=toremove[-1]
if(length(toremove)>0) {keep=keep[-toremove]}
}
Xw=Wd*X[index,keep]
if(!is.null(W)) {
    Xw=cbind(Xw,PhWy[index]*Wd)
	keep=c(keep,ncol(X)+1)
	}
if(remove_local_outlier){
cooks=cooks.distance(lm(YYw~Xw-1))
if(is.null(outv)) i_reject=which(cooks>4/length(YYw)) else i_reject=which(cooks>quantile(cooks,1-outv))
ii<-which(i_reject == which(index == i))
if(length(ii)>0) i_reject<-i_reject[-which(i_reject==which(index==i))]
if(length(i_reject)>0){
  Xw<-Xw[-i_reject,]
  YYw<-YYw[-i_reject]
  }
}
gwrm=fastlmLLT_C(as.matrix(Xw),as.matrix(YYw),SE)
Betav[i,keep]=gwrm$Betav
if(SE & !isgcv) {
  SEV[i,keep]=gwrm$se
  Zwi=solve(crossprod(Xw,Xw)) %*% t(Xw)
  tS=tS+(Xw[which(index==i),] %*% Zwi)[,which(index==i)]
}
if(betacor & !is.null(W)) if (Betav[i,ncol(X)+1]>1 | Betav[i,ncol(X)+1]< (-1)) Betav[i,keep]<-SARHS(as.matrix(Xw),as.matrix(YYw),SE,W[index,index],Betav[i,keep])
}
} else
{
registerDoParallel(ncore)
Betav<-foreach(i =1:n,.combine="rbind",.inorder=FALSE)  %dopar% {
betav<-rep(0,ifelse(is.null(W),m,m+1))
sev<-rep(0,ifelse(is.null(W),m,m+1))
if(is.null(Wh)) {
w.i<-GPKj(i-1,Y,X,S, H,kernels, type, minv = minv, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay) } else w.i<-Wh[i,]
index=which(w.i>0) #index=which((w.i[w.i>0]*length(w.i))>wimin)
Wd<-sqrt(w.i[index])
YYw<-Wd*Y[index]
keep=1:ncol(X)
if(m>1) {
toremove=unique(c(which(apply(X[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))),which((apply(Wd*X[index,],2,function(x) sum(abs(x)))/apply(XV[index,],2,function(x) sum(abs(x))))<xratiomin))) #toremove=which(apply(X[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
if(interceptv) toremove=toremove[-1]
if(length(toremove)>0) {keep=keep[-toremove]}
}
Xw=Wd*X[index,keep]
if(!is.null(W)) {
    Xw=cbind(Xw,PhWy[index]*Wd)
	keep=c(keep,ncol(X)+1)
	}
if(remove_local_outlier){
cooks=cooks.distance(lm(YYw~Xw-1))
if(is.null(outv)) i_reject=which(cooks>4/length(YYw)) else i_reject=which(cooks>quantile(cooks,1-outv))
ii<-which(i_reject == which(index == i))
if(length(ii)>0) i_reject<-i_reject[-which(i_reject==which(index==i))]
Xw<-Xw[-i_reject,]
YYw<-YYw[-i_reject]
}

gwrm=fastlmLLT_C(as.matrix(Xw),as.matrix(YYw),SE)
betav[keep]=gwrm$Betav
if(SE  & !isgcv) {
  sev[keep]=gwrm$se
  Zwi=solve(crossprod(Xw,Xw)) %*% t(Xw)
  tS=(Xw[which(index==i),] %*% Zwi)[,which(index==i)]
}
if(betacor & !is.null(W)) if (betav[ncol(X)+1]>1 | betav[ncol(X)+1]< (-1)) betav[keep]<-SARHS(as.matrix(Xw),as.matrix(YYw),SE,W[index,index],betav[keep])
if(!SE | isgcv) {betav} else {c(betav,sev,tS)}
}
if(!is.null(W)) m=m+1
if(SE & !isgcv) {SEV=Betav[,(m+1):(2*m)];tS=sum(Betav[,2*m+1])}
Betav=Betav[,1:m]
}
if(SE & !isgcv) list(Betav=Betav,SEV=SEV,edf=n-tS,tS=tS) else list(Betav=Betav,SEV=NULL,edf=NULL,tS=NULL)
}
