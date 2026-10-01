# =============================================================================
# Session caches of prep_w (weights) and MGWRSAR (distance maxima)
# =============================================================================
# The cached path must return exactly the values of the uncached path, for
# every Type, and must not reuse an entry after the distance matrix changed.

library(testthat)
library(mgwrsar)

prep_d <- mgwrsar:::prep_d
prep_w <- mgwrsar:::prep_w

set.seed(11)
n <- 120
co <- cbind(runif(n), runif(n), sample(1:365, n, replace = TRUE))

cases <- list(
  GD_gauss   = list(coords = co[, 1:2], kernels = "gauss", Type = "GD", H = 0.2, adaptive = FALSE, alpha = 1),
  GD_bisq_ad = list(coords = co[, 1:2], kernels = "bisq", Type = "GD", H = 25, adaptive = TRUE, alpha = 1),
  T_sym      = list(coords = co[, 3, drop = FALSE], kernels = "gauss_SYM_365", Type = "T", H = 30, adaptive = FALSE, alpha = 1),
  GDT_prod   = list(coords = co, kernels = c("gauss", "gauss_SYM_365"), Type = "GDT", H = c(0.2, 30), adaptive = c(FALSE, FALSE), alpha = 1),
  GDT_mixed  = list(coords = co, kernels = c("gauss", "gauss_SYM_365"), Type = "GDT", H = c(0.2, 30), adaptive = c(FALSE, FALSE), alpha = 0.5)
)

run_prep_w <- function(cs, dists, indexG) {
  prep_w(H = cs$H, kernels = cs$kernels, Type = cs$Type, adaptive = cs$adaptive,
         dists = dists, indexG = indexG, alpha = cs$alpha)$Wd
}

test_that("cached weights are identical to uncached weights (first call, stored, hit)", {
  old <- options(mgwrsar.weights_cache_mb = 256)
  on.exit(options(old), add = TRUE)
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    pd <- prep_d(coords = as.matrix(cs$coords), NN = n, TP = 1:n, extrapol = FALSE,
                 kernels = cs$kernels, Type = cs$Type)
    options(mgwrsar.weights_cache_mb = 0)
    ref <- run_prep_w(cs, pd$dists, pd$indexG)
    options(mgwrsar.weights_cache_mb = 256)
    mgwrsar:::.mgwrsar_clear_cache()
    expect_identical(run_prep_w(cs, pd$dists, pd$indexG), ref, info = paste(nm, "miss"))
    expect_identical(run_prep_w(cs, pd$dists, pd$indexG), ref, info = paste(nm, "stored"))
    expect_identical(run_prep_w(cs, pd$dists, pd$indexG), ref, info = paste(nm, "hit"))
    expect_true(abs(sum(ref) - nrow(ref)) < 1e-8, info = nm)
  }
})

test_that("a modified distance matrix is not served from the cache", {
  old <- options(mgwrsar.weights_cache_mb = 256)
  on.exit(options(old), add = TRUE)
  cs <- cases$GD_gauss
  pd <- prep_d(coords = as.matrix(cs$coords), NN = n, TP = 1:n, extrapol = FALSE,
               kernels = cs$kernels, Type = cs$Type)
  mgwrsar:::.mgwrsar_clear_cache()
  w1 <- run_prep_w(cs, pd$dists, pd$indexG)
  d2 <- pd$dists
  d2$dist_s[5, 7] <- d2$dist_s[5, 7] * 2          # copy-on-modify: new object
  options(mgwrsar.weights_cache_mb = 0)
  ref2 <- run_prep_w(cs, d2, pd$indexG)
  options(mgwrsar.weights_cache_mb = 256)
  expect_identical(run_prep_w(cs, d2, pd$indexG), ref2)
  expect_false(identical(ref2, w1))
  expect_identical(run_prep_w(cs, pd$dists, pd$indexG), w1)   # original still cached and valid
})

test_that("cached distance maxima follow the matrix they were computed on", {
  mgwrsar:::.mgwrsar_clear_cache()
  d <- matrix(c(0, 1, 2, 3), 2)
  expect_identical(mgwrsar:::.mgwrsar_dist_max(d), 3)
  expect_identical(mgwrsar:::.mgwrsar_dist_max(d), 3)
  d[2, 2] <- 10
  expect_identical(mgwrsar:::.mgwrsar_dist_max(d), 10)
  d_na <- matrix(c(NA, 1, 2, 4), 2)
  expect_identical(mgwrsar:::.mgwrsar_dist_max(d_na), 4)
})

test_that("the default cache budget scales with the matrix size", {
  old <- options(mgwrsar.weights_cache_mb = NULL)
  on.exit(options(old), add = TRUE)
  budget <- mgwrsar:::.mgwrsar_weights_budget
  expect_identical(budget(1000 * 1000), 512 * 2^20)          # 8 MB matrices: floor
  expect_identical(budget(4000 * 4000), 20 * 8 * 16e6)        # 128 MB matrices: 20 of them
  expect_identical(budget(10000 * 10000), 4096 * 2^20)        # 800 MB matrices: cap
  options(mgwrsar.weights_cache_mb = 64)
  expect_identical(budget(10000 * 10000), 64 * 2^20)          # explicit option wins
  options(mgwrsar.weights_cache_mb = 0)
  expect_identical(budget(1000 * 1000), 0)
})

test_that("column subsets are served from the cache and follow the matrix", {
  mgwrsar:::.mgwrsar_clear_cache()
  cols <- mgwrsar:::.mgwrsar_cols
  x <- matrix(runif(60), 6); storage.mode(x) <- "double"
  s1 <- cols(x, 4); s2 <- cols(x, 4)
  expect_identical(s1, x[, 1:4])
  expect_true(mgwrsar:::.mgwrsar_same_object(s1, s2))   # same object on the second call
  expect_identical(cols(x, 7), x[, 1:7])
  y <- x; y[2, 3] <- -1
  expect_identical(cols(y, 4), y[, 1:4])
  expect_false(identical(cols(y, 4), s1))
  xi <- matrix(1:60, 6)
  expect_identical(cols(xi, 5), xi[, 1:5])
})

test_that("scale truncation keeps every column that can carry a positive weight", {
  mgwrsar:::.mgwrsar_clear_cache()
  d <- t(apply(matrix(runif(40, 0, 1), 8), 1, sort)); d[, 1] <- 0     # kNN-like rows
  nc <- mgwrsar:::.mgwrsar_ncols_within
  for (h in c(0.05, 0.3, 0.6, 2)) {
    NN <- nc(d, h)
    expect_true(NN >= 3 && NN <= ncol(d))
    if (NN < ncol(d)) expect_true(all(d[, (NN + 1):ncol(d)] >= h))   # dropped columns: zero weights
  }
  expect_identical(nc(d, 2), ncol(d))
  sc <- mgwrsar:::.mgwrsar_scale_ncols
  ctrl <- list(adaptive = FALSE, dists = list(dist_s = d))
  expect_identical(sc(0.3, "bisq", ctrl, list()), nc(d, 0.3))
  expect_null(sc(0.3, "gauss", ctrl, list()))                                   # off by default
  expect_identical(sc(0.1, "gauss", ctrl, list(trunc_gauss = 3)), nc(d, 0.3))   # c * h
  expect_identical(sc(10, "bisq", list(adaptive = TRUE), list()), 12)           # adaptive compact: v + 2
  expect_null(sc(10, "gauss", list(adaptive = TRUE), list()))
})
