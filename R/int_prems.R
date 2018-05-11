#' int_prems
#' to be documented
#' @usage int_prems(X)
#' @param X  to be documented
#' @keywords internal
#' @return to be documented
int_prems <-
function(X){
col_int<-which(apply(X,2,function(x) sd(x)==1-mean(x)))
if(sum(col_int)!=0) cbind(X[,col_int],X[,-col_int]) else X
}
