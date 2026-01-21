#' reord_D
#'
#' Reorders the elements of a distance matrix row by row according to a corresponding
#' index matrix.
#'
#' @usage reord_D(M, idxG)
#' @param M A numeric matrix to be reordered.
#' @param idxG An integer matrix of the same number of rows as `M`, where
#' each row contains column indices indicating the new ordering for that row.
#' @return A reordered numeric matrix with the same number of rows as `M` and
#' columns equal to the maximum index in `idxG`.
#' @export
reord_D <- function(M, idxG) {
  .Call("_mgwrsar_knn_stable_sort", M, idxG, PACKAGE = "mgwrsar")
}
