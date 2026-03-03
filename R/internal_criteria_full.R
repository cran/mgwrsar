# Internal criteria helpers (FULL trace only)
# These helpers are intentionally non-exported. They centralize AICc/ENP/trace computations
# to reduce code drift across models (GWR/MGWR/TDS/ATDS; GD/GDT).
#
# IMPORTANT: In this audit track, we keep ONLY the "full" trace(S) mode.
# Any legacy block/band/sampling alternatives remain out of the criterion pipeline.

.compute_trace_full <- function(S) {
  sum(diag(S))
}

.compute_enp_from_trace <- function(trS) {
  as.numeric(trS)
}

# AICc from RSS and ENP (trace of S).
# Consistent with aicc_f(): n*log(RSS/n) + n*log(2*pi) + n*(n+enp)/(n-1-enp)
.compute_AICc_from_rss_enp <- function(n, rss, enp) {
  n <- as.numeric(n); rss <- as.numeric(rss); enp <- as.numeric(enp)

  if (!is.finite(n) || n <= 0) stop("n must be positive", call. = FALSE)
  if (!is.finite(rss) || rss <= 0) stop("rss must be positive and finite", call. = FALSE)
  if (!is.finite(enp) || enp < 0) stop("enp must be non-negative and finite", call. = FALSE)

  denom <- n - 1 - enp
  if (denom <= 0) return(Inf)
  n * log(rss / n) + n * log(2 * pi) + n * (n + enp) / denom
}
