# INTERNAL: thread control guardrails (OpenMP/BLAS)

.mgwrsar_is_check_context <- function() {
  val <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  val %in% c("true", "1", "yes")
}

.mgwrsar_normalize_nthreads <- function(nthreads = 1L) {
  nthreads <- suppressWarnings(as.integer(nthreads))
  if (!is.finite(nthreads) || nthreads < 1L) {
    nthreads <- 1L
  }
  if (.mgwrsar_is_check_context()) {
    nthreads <- 1L
  }
  as.integer(nthreads)
}

.mgwrsar_set_native_threads <- function(nthreads = 1L) {
  nthreads <- .mgwrsar_normalize_nthreads(nthreads)

  Sys.setenv(
    OMP_NUM_THREADS = as.character(nthreads),
    OPENBLAS_NUM_THREADS = as.character(nthreads),
    MKL_NUM_THREADS = as.character(nthreads),
    VECLIB_MAXIMUM_THREADS = as.character(nthreads),
    BLIS_NUM_THREADS = as.character(nthreads)
  )

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    # blas_set_num_threads costs ~0.3 ms and this helper runs before every
    # native routine (thousands of times in TDS_MGWR). The single-thread
    # request is skipped when this helper already applied it: the package's
    # own direct RhpcBLASctl calls only ever set 1 thread, so the cached state
    # cannot be stale for 1; any other count is always re-applied.
    if (!(nthreads == 1L && identical(.mgwrsar_thread_state$blas, 1L))) {
      try(RhpcBLASctl::blas_set_num_threads(nthreads), silent = TRUE)
      .mgwrsar_thread_state$blas <- nthreads
    }
    try(RhpcBLASctl::omp_set_num_threads(nthreads), silent = TRUE)
  }

  invisible(nthreads)
}

# Session cache: last BLAS thread count applied by .mgwrsar_set_native_threads,
# and number of physical cores (parallel::detectCores() spawns a system call).
.mgwrsar_thread_state <- new.env(parent = emptyenv())
.mgwrsar_thread_state$blas <- NA_integer_
.mgwrsar_thread_state$cores <- NA_integer_

.mgwrsar_physical_cores <- function() {
  if (is.na(.mgwrsar_thread_state$cores)) {
    dc <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (!is.finite(dc) || dc < 1L) dc <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (!is.finite(dc) || dc < 1L) dc <- 1L
    .mgwrsar_thread_state$cores <- as.integer(dc)
  }
  .mgwrsar_thread_state$cores
}
