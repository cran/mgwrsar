# INTERNAL: session caches of kernel weight matrices and distance maxima
#
# TDS_MGWR calls MGWRSAR thousands of times with the same distance matrices and
# a few dozen distinct bandwidths per axis, and prep_w used to recompute every
# kernel matrix from scratch. Entries are keyed by pointer identity of the
# distance matrix and keep a reference to it: while referenced, the matrix
# cannot be freed (so its address cannot be reused) and any R-level
# modification copies it (copy-on-modify), so a pointer match means the same,
# unchanged content. Cached results are therefore identical to recomputed ones.
#
# Memory: weight matrices requested at least twice are held, least recently
# used first out, up to getOption("mgwrsar.weights_cache_mb") megabytes (0
# disables the cache). When the option is unset, the budget holds 20 matrices
# of the current size, within [512 MB, 4 GB]: 512 MB up to n x NN ~ 3.4e6,
# 4 GB from n = NN ~ 5200 (at n = NN = 5000, 512 MB held only 2 matrices and
# tds_mgtwr took 579 s instead of 323 s).

.mgwrsar_cache <- new.env(parent = emptyenv())
.mgwrsar_cache$w <- list()
.mgwrsar_cache$seen <- list()
.mgwrsar_cache$max <- list()
.mgwrsar_cache$cols <- list()
.mgwrsar_cache$colmin <- list()

.mgwrsar_same_object <- function(a, b) {
  .Call("_mgwrsar_same_object", a, b, PACKAGE = "mgwrsar")
}

.mgwrsar_clear_cache <- function() {
  .mgwrsar_cache$w <- list()
  .mgwrsar_cache$seen <- list()
  .mgwrsar_cache$max <- list()
  .mgwrsar_cache$cols <- list()
  .mgwrsar_cache$colmin <- list()
  invisible(NULL)
}

# max(d, na.rm = TRUE) as a double, computed once per distance matrix
.mgwrsar_dist_max <- function(d) {
  for (e in .mgwrsar_cache$max) {
    if (.mgwrsar_same_object(e$d, d)) return(e$value)
  }
  v <- as.numeric(max(d, na.rm = TRUE))
  .mgwrsar_cache$max <- c(list(list(d = d, value = v)), utils::head(.mgwrsar_cache$max, 3L))
  v
}

# Cache budget in bytes for weight matrices of n_elem elements
.mgwrsar_weights_budget <- function(n_elem) {
  opt <- getOption("mgwrsar.weights_cache_mb")
  if (!is.null(opt)) return(opt * 2^20)
  min(4096 * 2^20, max(512 * 2^20, 20 * 8 * as.numeric(n_elem)))
}

# Weight matrix for distance matrix `d` and `key` (list of the other inputs),
# returned from the cache or computed by `compute()` and stored.
.mgwrsar_weights <- function(d, key, compute) {
  budget <- .mgwrsar_weights_budget(length(d))
  if (!is.numeric(budget) || length(budget) != 1L || is.na(budget) || budget <= 0) return(compute())
  ent <- .mgwrsar_cache$w
  for (i in seq_along(ent)) {
    e <- ent[[i]]
    if (identical(e$key, key) && .mgwrsar_same_object(e$d, d)) {
      if (i > 1L) .mgwrsar_cache$w <- c(ent[i], ent[-i])
      return(e$w)
    }
  }
  w <- compute()
  # Admission: a matrix is stored only when its key is requested a second
  # time. A bandwidth search visits each bandwidth once, and retaining its
  # matrices only costs memory and garbage collection; TDS_MGWR revisits the
  # same few dozen bandwidths thousands of times.
  seen <- .mgwrsar_cache$seen
  again <- FALSE
  for (i in seq_along(seen)) {
    if (identical(seen[[i]]$key, key) && .mgwrsar_same_object(seen[[i]]$d, d)) {
      again <- TRUE
      seen <- seen[-i]
      break
    }
  }
  if (!again) {
    .mgwrsar_cache$seen <- c(list(list(d = d, key = key)), utils::head(seen, 63L))
    return(w)
  }
  .mgwrsar_cache$seen <- seen
  bytes <- 8 * as.numeric(length(w))
  if (bytes <= budget) {
    ent <- c(list(list(d = d, key = key, w = w, bytes = bytes)), ent)
    ent <- ent[cumsum(vapply(ent, function(e) e$bytes, numeric(1))) <= budget]
    .mgwrsar_cache$w <- ent
  }
  w
}

# normW(Ws * Wt), fused natively when both are double matrices with the same
# attributes (same values as the R expression).
.mgwrsar_prod_normW <- function(Ws, Wt) {
  if (is.matrix(Ws) && is.matrix(Wt) && is.double(Ws) && is.double(Wt) &&
      identical(attributes(Ws), attributes(Wt))) {
    return(.Call("_mgwrsar_wprod_norm_cpp", Ws, Wt, PACKAGE = "mgwrsar"))
  }
  normW(Ws * Wt)
}

# x[, 1:NN] served from the cache when the same matrix was already cut at NN.
# update_opt() truncates the neighbour columns of dists and indexG at every
# candidate bandwidth (adaptive compact kernels); a fresh subset per call would
# defeat the weights cache, which is keyed by the identity of the distance
# matrix. Same budget rule as the weights (size of the full matrix).
.mgwrsar_cols <- function(x, NN) {
  budget <- .mgwrsar_weights_budget(length(x))
  if (!is.numeric(budget) || length(budget) != 1L || is.na(budget) || budget <= 0) return(x[, 1:NN])
  ent <- .mgwrsar_cache$cols
  for (i in seq_along(ent)) {
    e <- ent[[i]]
    if (e$NN == NN && .mgwrsar_same_object(e$x, x)) {
      if (i > 1L) .mgwrsar_cache$cols <- c(ent[i], ent[-i])
      return(e$sub)
    }
  }
  sub <- x[, 1:NN]
  bytes <- as.numeric(object.size(sub))
  if (bytes <= budget) {
    ent <- c(list(list(x = x, NN = NN, sub = sub, bytes = bytes)), ent)
    ent <- ent[cumsum(vapply(ent, function(e) e$bytes, numeric(1))) <= budget]
    .mgwrsar_cache$cols <- ent
  }
  sub
}

# Column minima of a distance matrix, computed once per matrix (cached by
# pointer identity, reference kept).
.mgwrsar_colmin <- function(d) {
  for (e in .mgwrsar_cache$colmin) {
    if (.mgwrsar_same_object(e$d, d)) return(e$value)
  }
  v <- suppressWarnings(apply(d, 2, min, na.rm = TRUE))
  .mgwrsar_cache$colmin <- c(list(list(d = d, value = v)), utils::head(.mgwrsar_cache$colmin, 3L))
  v
}

# Number of leading neighbour columns of `d` that can carry a distance < h
# (the last column whose minimum is below h; at least 3 so the subset stays
# a matrix). Every dropped column holds only distances >= h, i.e. zero weights
# for a compact kernel of bandwidth h, so truncating there is exact.
.mgwrsar_ncols_within <- function(d, h) {
  keep <- which(.mgwrsar_colmin(d) < h)
  max(3L, if (length(keep)) max(keep) else 0L)
}
