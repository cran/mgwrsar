#' reord_M
#'
#' Reorders the elements of a matrix row by row according to a corresponding
#' index matrix. Each row of the index matrix specifies the new column order
#' for the corresponding row of the input matrix.
#'
#' @param M A numeric matrix to be reordered.
#' @param idxG An integer matrix of the same number of rows as `M`, where
#' each row contains column indices indicating the new ordering for that row.
#' @return A reordered numeric matrix with the same number of rows as `M` and
#' columns equal to the maximum index in `idxG`.
#' @details
#' This function performs a scatter-like reordering: for each row `i`, the
#' values of `M[i, ]` are placed in the positions specified by `idxG[i, ]`.
#' It is particularly useful for spatial or neighborhood-based rearrangements
#' where each observation has its own local indexing of neighbors.
#' @examples
#' M <- matrix(1:9, nrow = 3, byrow = TRUE)
#' idxG <- matrix(c(3, 1, 2, 1, 2, 3, 2, 3, 1), nrow = 3, byrow = TRUE)
#' reord_M_R(M, idxG)
#' @export
reord_M_R <- function(M, idxG) {
  n <- nrow(idxG)
  k <- ncol(idxG)
  ncol_out <- max(idxG)

  # output initialisé avec fill
  out <- matrix(NA_real_, n, ncol_out)

  # Indices linéaires destinataires dans 'out'
  lin <- matrix(rep(1:n, each = k), nrow = n, byrow = TRUE)

  # Transformation en vecteur pour affectation rapide
  sel <- !is.na(idxG)
  out[cbind(lin[sel], idxG[sel])] <- M[sel]

  out
}
