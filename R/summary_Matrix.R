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
    attr(r, "header") <- sprintf("%d x %d sparse Matrix of class \"%s\", with %d entries",
        d[1], d[2], class(object), nnzero(object))
    class(r) <- c("sparseSummary", class(r))
    r
}
