#' normW
#' row normalization of dgCMatrix
#' @usage normW(x)
#' @param x  A dgCMatrix class matrix
#' @return A row normalized dgCMatrix
normW <- function(x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  rs <- Matrix::rowSums(x, na.rm = TRUE)
  rs[rs == 0] <- 1
  if (inherits(x, "sparseMatrix")) {
    return(Matrix::Diagonal(x = 1/rs) %*% x)
  } else {
    return(x / rs)
  }
}
