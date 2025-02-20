#' BP_pred_MGWRSAR
#' to be documented
#' @usage BP_pred_MGWRSAR(YS,X,W,e,beta_hat,lambda_hat,S,O,coords,
#' type='BPN',k=16,Wk=NULL)
#' @param YS to be documented
#' @param X to be documented
#' @param W to be documented
#' @param e to be documented
#' @param beta_hat to be documented
#' @param lambda_hat to be documented
#' @param S to be documented
#' @param O to be documented
#' @param coords to be documented
#' @param type to be documented
#' @param k to be documented
#' @param Wk to be documented
#' @noRd
#'
#' @return to be documented
BP_pred_MGWRSAR <-
function(YS,X,W,e,beta_hat,lambda_hat,S,O,coords,type='BPN',k=16,Wk=NULL){
###  YS = c(Y_insample)
###  X = rbind(X_insample,X_outsample)
###  e=residuals_insample
###  beta_hat,lambda_hat = coef_insample
###  S,O index of insample and outsample
## computation of Beta_hat on the whole sample
if(is.null(Wk)) Wk=(Beta_extropolation(coords=coords,k=k,O,type='S+OxS+O'))[-S,-O]
beta_hat=rbind(beta_hat,as.matrix(Wk %*% beta_hat))
if(length(lambda_hat)>1) lambda_hat=as.numeric(c(lambda_hat,as.numeric(Wk %*% lambda_hat)))
###
sigma2=mean(e^2)
J=S[which(Matrix::rowSums(W[S, O]) > 0)]
JO=unique(c(sort(J),sort(O)))
WJ<-W[JO,JO]
rs<-Matrix::rowSums(WJ)
sup0<-rs>0
WJ[sup0,]<-WJ[sup0,]/rs[sup0]
m=length(JO)

if(length(lambda_hat)>1) lambda_hat<-as.numeric(lambda_hat[JO])
### creation des index sur J et O
J2=1:length(J)
O2=(length(J)+1):m
## Compute Q
Q=(Diagonal(m,1)-lambda_hat*(WJ+Matrix::t(WJ))+lambda_hat^2*(Matrix::t(WJ)%*%WJ))/sigma2
## Compute XB
XB=Matrix::rowSums(X[JO,]* beta_hat[JO,])
## Compute TC
if(m<2000) iWJ=solve(Diagonal(m,1)-lambda_hat*WJ) else iWJ=ApproxiW(WJ,lambda_hat,8)
YTC=iWJ %*% XB
## compute BP
YBPN=YTC[O2]-solve(Q[O2,O2])%*%Q[O2,J2]%*% (YS[J]-YTC[J2])
if(type=='BPN') res=YBPN else res=YTC[O2]
res
}
