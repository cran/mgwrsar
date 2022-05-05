#' comb
#' combine method for foreach
#' @usage comb(...)
#' @return combined list
#' @noRd
comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}
