# INTERNAL: package-internal random numbers without touching the user's generator
#
# MGWRSAR(), the bandwidth searches, multiscale_gwr(), TDS_MGWR(),
# simu_multiscale() and the bootstrap test seed an L'Ecuyer-CMRG generator for
# their own draws (jitter of duplicated coordinates, sub-sampling of distances,
# bootstrap resamples, parallel streams). Until 1.4 this changed the generator
# of the user's session as a side effect: every set.seed() after a fit drew
# from L'Ecuyer-CMRG instead of the session's generator, so simulations
# around the package were not reproducible. The internal seeding is kept, so
# the results are unchanged; the user's generator (kind and state) is saved
# at entry and restored when the function returns.
#
# Usage, at the top of a function that seeds:
#   rng_state <- .mgwrsar_rng_save(); on.exit(.mgwrsar_rng_restore(rng_state), add = TRUE)

.mgwrsar_rng_save <- function() {
  list(kind = RNGkind(),
       seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
         get(".Random.seed", envir = globalenv()) else NULL)
}

.mgwrsar_rng_restore <- function(state) {
  suppressWarnings(do.call(RNGkind, as.list(state$kind)))
  if (is.null(state$seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", state$seed, envir = globalenv())
  }
  invisible(NULL)
}
