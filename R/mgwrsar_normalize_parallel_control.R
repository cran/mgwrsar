# INTERNAL: normalize parallel controls
.mgwrsar_normalize_parallel_control <- function(control, context = c("MGWRSAR", "search_bandwidths")) {
  context <- match.arg(context)

  if (is.null(control)) control <- list()

  # During checks/install contexts that limit cores, stay strictly single-core.
  limit_cores <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  if (limit_cores %in% c("true", "1", "yes")) {
    control$ncore <- 1L
    if ("doMC" %in% names(control)) control$doMC <- FALSE
    return(control)
  }

  has_ncore <- "ncore" %in% names(control)
  has_doMC <- "doMC" %in% names(control)

  # Backward compatibility: if doMC=TRUE and ncore missing, use available cores - 1
  if (has_doMC && isTRUE(control$doMC) && !has_ncore) {
    dc <- .mgwrsar_physical_cores()
    control$ncore <- max(1L, as.integer(dc - 1L))
    has_ncore <- TRUE
  }

  ncore_val <- if (has_ncore) suppressWarnings(as.integer(control$ncore)) else 1L
  if (!is.finite(ncore_val) || ncore_val < 1L) ncore_val <- 1L

  max_core <- .mgwrsar_physical_cores()

  if (ncore_val > max_core) {
    warning(
      sprintf("Requested ncore=%d exceeds available cores (%d). Using ncore=%d.", ncore_val, max_core, max_core),
      call. = FALSE
    )
    ncore_val <- max_core
  }

  control$ncore <- .mgwrsar_normalize_nthreads(ncore_val)
  control
}
