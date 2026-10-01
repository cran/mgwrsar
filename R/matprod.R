# INTERNAL: dense matrix product of the hat-matrix recursions (S %*% R_k, ...)
#
# Replaces SMUT::eigenMapMatMult by the same Eigen product compiled in the
# package. With an optimised BLAS, R's own %*% is much faster than Eigen's
# single-threaded product, while the reference BLAS shipped with R is much
# slower (n = 4000, Apple M-series: 0.28 s with Accelerate, 3.8 s with Eigen,
# 21 s with the reference BLAS). The backend therefore follows the BLAS that R
# is linked to: %*% when it is a known optimised library, Eigen otherwise.
# options(mgwrsar.matprod = "eigen") or "blas" forces one of them.

.mgwrsar_matprod_env <- new.env(parent = emptyenv())

.mgwrsar_blas_is_optimised <- function() {
  if (is.null(.mgwrsar_matprod_env$optimised)) {
    blas <- tryCatch(extSoftVersion()[["BLAS"]], error = function(e) "")
    .mgwrsar_matprod_env$optimised <-
      grepl("veclib|accelerate|openblas|mkl|blis|atlas|flexiblas|armpl", blas,
            ignore.case = TRUE)
  }
  .mgwrsar_matprod_env$optimised
}

.mgwrsar_matprod <- function(A, B) {
  backend <- getOption("mgwrsar.matprod", "auto")
  if (identical(backend, "blas") ||
      (!identical(backend, "eigen") && .mgwrsar_blas_is_optimised())) {
    C <- A %*% B
    dimnames(C) <- NULL
    return(C)
  }
  .Call("_mgwrsar_matprod_C", A, B, PACKAGE = "mgwrsar")
}
