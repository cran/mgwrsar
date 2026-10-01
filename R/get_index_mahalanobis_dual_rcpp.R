#' to be documented
#' @usage get_index_mahalanobis_dual_rcpp(DS, DTM)
#' @keywords internal
#' @return to be documented
#' @noRd
get_index_mahalanobis_dual_rcpp <- function(DS, DTM, nthreads = 1L) {
  .mgwrsar_set_native_threads(nthreads)
  .Call("_mgwrsar_get_index_mahalanobis_dual_rcpp", DS, DTM, PACKAGE = "mgwrsar")
}
