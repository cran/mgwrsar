# INTERNAL: normalize deprecated parallel controls
.mgwrsar_normalize_parallel_control <- function(control, context = c("MGWRSAR", "search_bandwidths")) {
  context <- match.arg(context)

  if (is.null(control)) control <- list()

  has_ncore <- "ncore" %in% names(control)

  ncore_val <- if (has_ncore) suppressWarnings(as.integer(control$ncore)) else 1L
  if (is.na(ncore_val) || ncore_val < 1L) ncore_val <- 1L

  # Condition that indicates user is asking for parallel estimation (deprecated)
  asked_parallel_estimation <- (has_ncore && ncore_val > 1L)

  if (asked_parallel_estimation) {
    # Warn only in MGWRSAR (estimation), not in search (where parallel still exists)
    if (identical(context, "MGWRSAR")) {
      msg <- paste0(
        "Parallel estimation in `MGWRSAR()` is no longer supported (arguments `control$doMC`/`control$ncore`). ",
        "Parallelism is now used only for bandwidth search via `search_bandwidths()` ",
        "for GWR/GTWR-type models. Proceeding with `ncore = 1`."
      )

      # avoid spamming: warn once per session
      if (!isTRUE(getOption("mgwrsar.warned_parallel_estimation", FALSE))) {
        options(mgwrsar.warned_parallel_estimation = TRUE)
        warning(msg, call. = FALSE)
      }
    }
  }

  # Always force-off parallel estimation for safety
  control$ncore <- 1L

  control
}
