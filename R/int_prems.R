#' Reorder Intercept Column First
#'
#' Identifies the intercept column (constant column) in a matrix and moves it
#' to the first position.
#'
#' @usage int_prems(X)
#' @param X A numeric matrix, typically the design matrix.
#' @keywords internal
#' @return The matrix \code{X} with the intercept column moved to the first position.
#'   If no intercept is found, returns \code{X} unchanged.
#' @export
int_prems <-
function(X){
col_int<-which(apply(X,2,function(x) sd(x)==1-mean(x)))
if(sum(col_int)!=0) cbind(X[,col_int],X[,-col_int]) else X
}
