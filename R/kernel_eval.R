#' Evaluate a kernel on a distance matrix
#'
#' Computes the kernel weights natively for bisq, gauss, epane and triangle
#' (fixed and \code{_adapt_sorted}), with the same values as the R kernels,
#' and falls back to the R kernel function for any other kernel or input the
#' native path does not reproduce exactly.
#'
#' @usage kernel_eval(kernel, d, h, n_norm = 0L)
#' @param kernel Name of the kernel function (e.g. \code{"bisq_adapt_sorted"}).
#' @param d Matrix of distances.
#' @param h Bandwidth (number of neighbours for adaptive kernels).
#' @param n_norm Number of \code{normW()} row normalisations applied to the
#'   kernel matrix (done in place on the native path, one allocation).
#' @return A matrix of kernel weights with the dimensions of \code{d}.
#' @keywords internal
#' @noRd
kernel_eval <- function(kernel, d, h, n_norm = 0L) {
  if (is.character(kernel) && length(kernel) == 1L &&
      is.matrix(d) && is.double(d) && is.numeric(h)) {
    w <- .Call("_mgwrsar_kernel_w_cpp", d, as.double(h), kernel, as.integer(n_norm),
               PACKAGE = "mgwrsar")
    if (!is.null(w)) return(w)
  }
  w <- do.call(kernel, args = list(d, h))
  for (r in seq_len(n_norm)) w <- normW(w)
  w
}
