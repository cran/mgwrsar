# INTERNAL: dense matrix product of the hat-matrix recursions (S %*% R_k, ...)
#
# Replaces SMUT::eigenMapMatMult by the same Eigen product compiled in the
# package. With an optimised BLAS, R's own %*% is much faster than Eigen's
# single-threaded product, while the reference BLAS shipped with R is much
# slower (n = 4000, Apple M-series: 0.28 s with Accelerate, 3.8 s with Eigen,
# 21 s with the reference BLAS). The backend therefore follows the BLAS that R
# is linked to: %*% when it is a known optimised library, Eigen otherwise.
# options(mgwrsar.matprod = "eigen") or "blas" forces one of them.
#
# The %*% result is verified (Freivalds check, O(n^2)) and recomputed when
# wrong. On Apple Silicon (M5 Max, macOS 26.6, Accelerate) the first large
# dgemm run by the R process within a few milliseconds after
# parallel::mclapply() returns is sporadically wrong (relative errors 1e-4 to
# 0.4 on the whole product), whatever VECLIB_MAXIMUM_THREADS. TDS_MGWR runs
# the hat-matrix recursion right after its parallel candidate search, so Shat,
# R_k, AICc and the standard errors varied from run to run with ncore > 1. A
# correct product passes the check by ten orders of magnitude (relative
# discrepancy < 1e-13 against >= 1e-4 for the corrupted ones).

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
    # A second %*% runs outside the faulty window; Eigen as a last resort
    for (attempt in 1:2) {
      C <- A %*% B
      dimnames(C) <- NULL
      if (.mgwrsar_matprod_ok(A, B, C)) return(C)
    }
  }
  .Call("_mgwrsar_matprod_C", A, B, PACKAGE = "mgwrsar")
}

# Freivalds check of C = A %*% B with a fixed positive vector (no RNG draw).
# Products with non-finite values cannot be checked and are accepted.
.mgwrsar_matprod_ok <- function(A, B, C) {
  if (length(C) == 0L) return(TRUE)
  v <- 1 + (seq_len(ncol(B)) %% 7) / 8
  r <- A %*% (B %*% v)
  d <- C %*% v - r
  if (!all(is.finite(r)) || !all(is.finite(d))) return(TRUE)
  max(abs(d)) <= 1e-8 * max(abs(r))
}
