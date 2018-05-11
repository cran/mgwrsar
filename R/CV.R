#' CV
#' to be documented
#' @usage CV(model)
#' @param model to be documented
#' @keywords internal
#'
#' @return to be documented
CV <-
function(model){
cat('\n ###### compute CV for ',model$Model,'model ... could be long\n')
formula=model$formula
data=model$data
Y=model$Y
X=model$X
W=model$W
Model=model$Model
n=nrow(data)
res=c()
for(i in 1:n){
if(Model=='OLS') {res[i]<-Y[i]-predict(lm(formula,data[-i,]),newdata=data[i,])} else {
model$Betac <- mod(as.matrix(Y[-i]),as.matrix(X[-i,]),W[-i,-i],as.matrix(X[-i,]),as.matrix(Y[-i]),rep(1,n-1),'L0',FALSE, FALSE)
iW<-solve(Diagonal(n,1)-model$Betac[length(model$Betac)]*W)
res[i]<-Y[i]-sum((iW %*% (X %*% model$Betac[-length(model$Betac)]))[i,]) ### A VERIFIER
}
}
list(CV=sum(res^2),residuals_i=res)
}
