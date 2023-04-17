#' summary_Matrix
#' to be documented
#' @usage summary_Matrix(object, ...)
#' @param object  to be documented
#' @param ...  to be documented
#' @return to be documented
summary_Matrix <-
function (object, ...)
{
    d <- dim(object)
    T <- as(object, "TsparseMatrix")
    r <- data.frame(i = T@i + 1L, j = T@j + 1L, x = T@x)
    attr(r, "header") <- paste0(d[1]," x ",d[2], "sparse Matrix of class '",class(object),"' with",nnzero(object))
    class(r) <- c("sparseSummary", class(r))
    r
}
