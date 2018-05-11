#' MGWR
#' to be documented
#' @usage MGWR(Y,XC,XV,S,H, kernels, type = "GD",model='MGWR', minv = 1,
#' maxknn = 500, NmaxDist = 6000,SE=FALSE, isgcv, TIME, decay,
#' interceptv=TRUE,W=NULL,betacor=FALSE,remove_local_outlier=FALSE,
#' outv=0,doMC=FALSE,ncore=1,Wh=NULL,xratiomin=10e-10)
#' @param Y  A vector
#' @param XC  A matrix with covariates with stationnary parameters
#' @param XV   A matrix with covariates with spatially varying parameters
#' @param S  A matrix with variables used in kernel
#' @param H  A vector of bandwidths
#' @param kernels  A vector of kernel types
#' @param type  Type of Genelarized kernel product ('GD' only spatial,'GDC' spatial + a categorical variable,'GDX' spatial + a continuous variable,'GDT' spatial + a time index, and other combination 'GDXXC','GDTX',...)
#' @param model  A mgwrsar model type (see MGWRSAR)
#' @param minv  Minimum number of non null weight
#' @param NmaxDist  Maximum number of observation for computing dense weight matrix
#' @param maxknn  If n >NmaxDist how many column with dense weight matrix (max number of neighbours)
#' @param SE  If standard error are computed
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
MGWR <-
function(Y,XC,XV,S,H, kernels, type = "GD",model='MGWR', minv = 1, maxknn = 500, NmaxDist = 6000,SE=FALSE, isgcv, TIME, decay,interceptv=TRUE,W=NULL,betacor=FALSE,remove_local_outlier=FALSE,outv=0,doMC=FALSE,ncore=1,Wh=NULL,xratiomin=10e-10){
se=NULL
sev=NULL
ptm<-proc.time()
if(!is.null(XC)) XCG<-XC<-as.matrix(XC)
if(!is.null(XV)) XVG<-XV<-as.matrix(XV)
X=cbind(XC,XV)
WY=as.matrix(W%*%Y)

if(model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv')){W2=W
XVG<-cbind(XVG,WY)
} else W2=NULL
n<- NROW(Y)
m<-ncol(XV)
K<-ncol(XC)
SY<- matrix(0,nrow=n, ncol=1)

## model 'MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv'
if(model %in% c('MGWRSAR_0_kc_kv','MGWRSAR_0_0_kv')){
PhWy=PhWY_C(as.matrix(Y),as.matrix(X),W,rep(1,n))
if(model =='MGWRSAR_0_kc_kv') XC=cbind(XC,PhWy) else XC=PhWy
if(model =='MGWRSAR_0_kc_kv') XCG=cbind(XC,WY) else XCG=WY
}
XX<- XC

## step1
if(doMC==FALSE){
for (i in 1:n){
if(is.null(Wh)) {
w.i<-GPKj(i-1,Y,XV,S, H, kernels, type, minv = minv, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay) } else w.i<-Wh[i,]
index=which(w.i>0) #index=which((w.i[w.i>0]*length(w.i))>wimin)
wi<-sqrt(w.i)
Wd<-wi[index]
Yw<-Wd*Y[index]
keep=1:ncol(XV)
if(m>1) {
toremove=unique(c(which(apply(XV[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))),which((apply(Wd*XV[index,],2,function(x) sum(abs(x)))/apply(XV[index,],2,function(x) sum(abs(x))))<xratiomin))) #toremove=which(apply(XV[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
if(interceptv) toremove=toremove[-1]
if(length(toremove)>0) {
    	keep=keep[-toremove]
    	}
}
Xw=as.matrix(Wd*XV[index,keep])
if(model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv')){
    PhWy=PhWY_C(as.matrix(Y[index]),as.matrix(X[index,]*Wd),W[index,index],rep(1,length(index)))*Wd
Xw=cbind(Xw,PhWy)
XVi<-c(XV[i,keep],WY[i])
} else XVi<-XV[i,keep]

### todo : remove outlier for first stage, this first try doesn't work
# if(remove_local_outlier){
#   cooks=cooks.distance(lm(Y[index]~Xw-1))
#   if(is.null(outv)) i_reject=which(cooks>4/length(index)) else i_reject=which(cooks>quantile(cooks,1-outv))
#   ii<-which(i_reject == which(index == i))
#   if(length(ii)>0) i_reject<-i_reject[-which(i_reject==which(index==i))]
#   Xw<-Xw[-i_reject,]
#   matB<-QRcpp2_C(Xw[-i_reject,],as.matrix(Wd[-i_reject]*Y[index[-i_reject]]),as.matrix(Wd[-i_reject]*XC[index[-i_reject],]))
# } else
matB<-QRcpp2_C(Xw,as.matrix(Wd*Y[index]),as.matrix(Wd*XC[index,]))


XX[i,]<- XX[i,]-XVi %*% matB$XCw ## keep 2 option2 : =  XX[i,]-XV[i,keep] %*% qr.solve(Xw,Wd*XC[index,])
													      ##			 ou   = -XV[i,keep] %*% qr.solve(Xw,Wd*XC[index,])
SY[i]<- Y[i]-XVi %*% matB$SY  ## keep 2 option2 : =  Y[i]-XV[i,keep] %*% qr.solve(Xw,Wd*Y[index,])
}
} else {
registerDoParallel(ncore)
XXS<-foreach(i =1:n,.combine="rbind",.inorder=FALSE)  %dopar% {
if(is.null(Wh)) {
w.i<-GPKj(i-1,Y,XV,S, H, kernels, type, minv = minv, maxknn = 500, NmaxDist = 6000, isgcv, TIME, decay) } else w.i<-Wh[i,]
index=which(w.i>0) # index=which((w.i[w.i>0]*length(w.i))>wimin)
wi<-sqrt(w.i)
Wd<-wi[index]
Yw<-Wd*Y[index]
keep=1:ncol(XV)
if(m>1) {
toremove=unique(c(which(apply(XV[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))),which((apply(Wd*XV[index,],2,function(x) sum(abs(x)))/apply(XV[index,],2,function(x) sum(abs(x))))<xratiomin))) #toremove=which(apply(XV[index,], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
if(interceptv) toremove=toremove[-1]
if(length(toremove)>0) {
    	keep=keep[-toremove]
    	}
}
Xw=as.matrix(Wd*XV[index,keep])
if(model %in% c('MGWRSAR_1_0_kv','MGWRSAR_1_kc_kv')){
    PhWy=PhWY_C(as.matrix(Y[index]),as.matrix(X[index,]*Wd),W[index,index],rep(1,length(index)))*Wd
Xw=cbind(Xw,PhWy)
XVi<-c(XV[i,keep],WY[i])
} else XVi<-XV[i,keep]

### todo : remove outlier for first stage, this first try doesn't work
# if(remove_local_outlier){
#   cooks=cooks.distance(lm(Y[index]~Xw-1))
#   if(is.null(outv)) i_reject=which(cooks>4/length(index)) else i_reject=which(cooks>quantile(cooks,1-outv))
#   ii<-which(i_reject == which(index == i))
#   if(length(ii)>0) i_reject<-i_reject[-which(i_reject==which(index==i))]
#   Xw<-Xw[-i_reject,]
#   matB<-QRcpp2_C(Xw[-i_reject,],as.matrix(Wd[-i_reject]*Y[index[-i_reject]]),as.matrix(Wd[-i_reject]*XC[index[-i_reject],]))
# } else

  matB<-QRcpp2_C(Xw,as.matrix(Wd*Y[index]),as.matrix(Wd*XC[index,]))

XXt<- XX[i,]-XVi %*% matB$XCw ## keep 2 option2 : =  XX[i,]-XV[i,keep] %*% qr.solve(Xw,Wd*XC[index,])
													      ##			 ou   = -XV[i,keep] %*% qr.solve(Xw,Wd*XC[index,])
SYt<- Y[i]-XVi %*% matB$SY  ## keep 2 option2 : =  Y[i]-XV[i,keep] %*% qr.solve(Xw,Wd*Y[index,])
cbind(XXt,SYt)
}
XX=XXS[,1:(ncol(XXS)-1)]
SY=XXS[,ncol(XXS)]
}
if(!SE){
Betac<-QRcpp_C(as.matrix(XX),as.matrix(SY)) ### return1
} else {
modelsar<-fastlmLLT_C(as.matrix(XX),as.matrix(SY),SE)
Betac=modelsar$Betav
se=modelsar$se
}
YY<-Y-XC %*%as.matrix(Betac) ### return2

## step2

modelGWR<-GWR(YY,XV,X,Y,S,H, kernels, type, minv, maxknn, NmaxDist,SE, isgcv, TIME,decay,interceptv,W=W2,betacor,remove_local_outlier=remove_local_outlier,
outv=outv,doMC=doMC,ncore=ncore,Wh=Wh)
Betav=modelGWR$Betav
edf=modelGWR$edf-ncol(XC)
tS=modelGWR$tS
if(SE) sev=modelGWR$SEV
#Betav=(GWRSAR_C(as.matrix(YY),as.matrix(XV),W,as.matrix(coord),as.matrix(X),as.matrix(Y),71,FALSE,'bisq_knn','D', 100, 50000, FALSE,0,model,'2SLS','L5',FALSE))$Betav
diff<-proc.time()-ptm
#cat("\n Duree",diff[3], " secondes")
z <- list(Betav=Betav,Betac=Betac,se=se,SEV=sev,XV=XVG,XC=XC,edf=edf,tS=tS)
invisible(z)
}
