# =============================================================================
# Tests for the native local variance kernel (src/local_var.cpp) used by
# get_HKmin() in TDS_MGWR() when control_tds$check_pairs = TRUE
# =============================================================================
# Package: mgwrsar
#
# What is checked:
# - local_var_cpp(Wd, indexG, X, nthreads) equals the R reference that gathers
#   X through the neighbour index (prep_w() weights are n x NN in neighbour
#   order), for GD and GDT weights, with NN = n and NN < n, to 1e-10
# - the thread count does not change the result (sequential without OpenMP)
# - Wd with fewer columns than indexG (truncated kernel) is accepted and only
#   its columns are read; inconsistent inputs are rejected
# =============================================================================

library(testthat)
library(mgwrsar)

lv_n <- 100L
set.seed(11)
lv_co <- cbind(runif(lv_n), runif(lv_n))
lv_tm <- runif(lv_n, 1, 50)
lv_X <- cbind(a = rnorm(lv_n), b = lv_co[, 1] + 0.1 * rnorm(lv_n), c = as.numeric(lv_co[, 2] > 0.5))

# R reference: for each site i, weighted moments over its neighbours
lv_ref <- function(Wd, indexG, X) {
  m <- ncol(Wd)
  t(sapply(seq_len(nrow(Wd)), function(i) {
    xi <- X[indexG[i, seq_len(m)], , drop = FALSE]
    colSums(Wd[i, ] * xi^2) - colSums(Wd[i, ] * xi)^2
  }))
}
lv_cpp <- function(Wd, indexG, X, nthreads = 1L)
  .Call("_mgwrsar_local_var_cpp", Wd, indexG, X, as.integer(nthreads), PACKAGE = "mgwrsar")

test_that("local_var_cpp matches the R gather reference (GD, NN = n and NN < n)", {
  for (NN in c(lv_n, 30L)) {
    G <- mgwrsar:::prep_d(coords = lv_co, NN = NN, TP = 1:lv_n, kernels = "bisq", Type = "GD")
    for (h in c(8, 25)) {
      Wd <- mgwrsar:::prep_w(H = h, kernels = "bisq", Type = "GD", adaptive = TRUE,
                             dists = G$dists, indexG = G$indexG)$Wd
      expect_equal(dim(Wd), c(lv_n, NN))
      ref <- lv_ref(Wd, G$indexG, lv_X)
      expect_equal(unname(lv_cpp(Wd, G$indexG, lv_X)), unname(ref), tolerance = 1e-10)
      expect_equal(lv_cpp(Wd, G$indexG, lv_X, 2L), lv_cpp(Wd, G$indexG, lv_X, 1L))
    }
  }
})

test_that("local_var_cpp matches the R gather reference (GDT, NN = n and NN < n)", {
  for (NN in c(lv_n, 30L)) {
    G <- mgwrsar:::prep_d(coords = cbind(lv_co, lv_tm), NN = NN, TP = 1:lv_n,
                          kernels = c("bisq", "gauss"), Type = "GDT")
    Wd <- mgwrsar:::prep_w(H = c(20, 10), kernels = c("bisq", "gauss"), Type = "GDT",
                           adaptive = c(TRUE, FALSE), dists = G$dists, indexG = G$indexG)$Wd
    expect_equal(dim(Wd), c(lv_n, NN))
    ref <- lv_ref(Wd, G$indexG, lv_X)
    expect_equal(unname(lv_cpp(Wd, G$indexG, lv_X)), unname(ref), tolerance = 1e-10)
  }
})

test_that("local_var_cpp reads only ncol(Wd) columns of indexG and rejects bad inputs", {
  G <- mgwrsar:::prep_d(coords = lv_co, NN = lv_n, TP = 1:lv_n, kernels = "bisq", Type = "GD")
  Wd <- mgwrsar:::prep_w(H = 8, kernels = "bisq", Type = "GD", adaptive = TRUE,
                         dists = G$dists, indexG = G$indexG)$Wd
  Wt <- Wd[, 1:20]   # truncated weights, full index
  expect_equal(unname(lv_cpp(Wt, G$indexG, lv_X)), unname(lv_ref(Wt, G$indexG, lv_X)), tolerance = 1e-10)
  expect_error(lv_cpp(Wd, G$indexG[, 1:10], lv_X), "at least ncol")
  expect_error(lv_cpp(Wd, G$indexG[1:50, ], lv_X), "nrow")
  bad <- G$indexG; bad[1, 1] <- lv_n + 1L
  expect_error(lv_cpp(Wd, bad, lv_X), "1..nrow")
})
