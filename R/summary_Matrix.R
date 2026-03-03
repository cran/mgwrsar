#' Summary of a Sparse Matrix
#'
#' Produces a summary of a sparse matrix object showing its non-zero elements
#' in triplet form (row, column, value).
#'
#' @usage summary_Matrix(object, ...)
#' @param object A sparse matrix object (inheriting from \code{Matrix}).
#' @param ... Additional arguments (currently unused).
#' @return A data frame with columns \code{i}, \code{j}, \code{x} representing the
#'   row indices, column indices, and values of non-zero elements.
#' @export
summary_Matrix <-
function (object, ...)
{
    d <- dim(object)
    Tsp <- as(object, "TsparseMatrix")
    r <- data.frame(i = Tsp@i + 1L, j = Tsp@j + 1L, x = Tsp@x)
    attr(r, "header") <- paste0(d[1]," x ",d[2], "sparse Matrix of class '",class(object),"' with",nnzero(object))
    class(r) <- c("sparseSummary", class(r))
    r
}
