#' normW
#' row normalization of dgCMatrix
#' @usage normW(W)
#' @param W  A dgCMatrix class matrix
#' @return A row normalized dgCMatrix
normW <-function (W)
{
    if (is(W, "matrix"))
        W <- Matrix(W)
    W <- drop0(W)
    ret <- summary_Matrix(as(W, "dgCMatrix"))
    ti <- tapply(ret$x, ret$i, function(x) sum(x, na.rm = TRUE))
    ret$x <- as.numeric(ret$x/ti[match(ret$i, as.numeric(names(ti)))])
    W <- sparseMatrix(i = ret$i, j = ret$j, x = ret$x, dims = dim(W))
    W[is.na(W)] <- 0
    W
}
