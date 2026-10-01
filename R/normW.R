#' Row Normalization of Sparse Matrix
#'
#' Row-normalizes a sparse matrix (dgCMatrix) or dense matrix so that each
#' row sums to 1. Rows with zero sums are left unchanged.
#'
#' @usage normW(x)
#' @param x A matrix or dgCMatrix class sparse matrix.
#' @return A row-normalized matrix of the same class as the input.
#' @export
normW <- function(x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  # dense double matrix: native path, same values as x / rowSums(x, na.rm = TRUE)
  if (is.matrix(x) && is.double(x)) {
    return(.Call("_mgwrsar_normw_dense_cpp", x, PACKAGE = "mgwrsar"))
  }
  rs <- Matrix::rowSums(x, na.rm = TRUE)
  rs[rs == 0] <- 1
  if (inherits(x, "sparseMatrix")) {
    return(Matrix::Diagonal(x = 1/rs) %*% x)
  } else {
    return(x / rs)
  }
}
