# =============================================================================
# Tests for prep_w() - Weight computation from distances/indexG
# =============================================================================
# Package: mgwrsar v1.3.2
# Function: prep_w()
#
# Preconditions:
# - prep_d() already validated (distances + indexG are correct)
# - Here we test only the transformation dists -> Wd
# =============================================================================

library(testthat)
library(mgwrsar)
library(stringr)

# -------------------------------------------------------------------------
# Assumes you already created:
# test_data_full <- setup_test_data_full()
# test_data      <- get_subset_30obs(test_data_full)
# and your corrected helpers exist in the test file (or sourced)
# -------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Helpers for prep_w tests
# -----------------------------------------------------------------------------

check_Wd_basic <- function(Wd, tol = 1e-12) {
  expect_true(is.matrix(Wd))
  expect_true(all(is.finite(Wd)))
  expect_true(all(Wd >= -tol))
  rs <- rowSums(Wd)
  expect_true(all(abs(rs - 1) < 1e-8))
}

# For gaussian-like kernels, closer distances should get >= weight than farther
check_monotone_decreasing_proxy <- function(dist_mat, Wd, tol = 1e-12) {
  expect_equal(dim(dist_mat), dim(Wd))
  i_min <- apply(dist_mat, 1, which.min)
  i_max <- apply(dist_mat, 1, which.max)
  w_min <- Wd[cbind(seq_len(nrow(Wd)), i_min)]
  w_max <- Wd[cbind(seq_len(nrow(Wd)), i_max)]
  expect_true(all(w_min + tol >= w_max))
}

# "Oracle" replication of prep_w logic for GD/T/GDT (for non-adaptive gauss)
oracle_prep_w_GD <- function(dists, H) {
  normW(do.call("gauss", list(dists[["dist_s"]], H)))
}

oracle_prep_w_T  <- function(dists, kernel_name, H) {
  # kernel_name is already extracted e.g. "gauss" or "gauss_adapt_sorted"
  normW(do.call(kernel_name, list(dists[["dist_t"]], H)))
}

oracle_prep_w_GDT <- function(dists, H, alpha = 1,
                              kernel_s = "gauss", kernel_t = "gauss") {
  Ws <- normW(do.call(kernel_s, list(dists[["dist_s"]], H[1])))
  Wt <- normW(do.call(kernel_t, list(dists[["dist_t"]], H[2])))
  
  if (alpha == 1) {
    W <- Ws * Wt
  } else if (alpha == 0) {
    W <- Ws + Wt
  } else {
    Term_Inter <- Ws * Wt
    W <- (Ws + Wt) * (1 - alpha) + Term_Inter * alpha
  }
  normW(W)
}

test_data <- setup_test_data_full(
  n          = 30,
  lambda     = NULL,
  config_beta= 'default',
  config_snr = 0.9
)
# =============================================================================
# GROUP 1: GD (spatial)
# =============================================================================

test_that("prep_w: GD gauss, NN=n, weights are valid and monotone proxy holds", {
  
  
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 1:n
  
  pd <- prep_d(coords = data$coords, NN = n, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = "gauss", Type = "GD")
  
  pw <- prep_w(H = c(n), kernels = "gauss", Type = "GD",
               adaptive = FALSE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  check_Wd_basic(pw$Wd)
  check_monotone_decreasing_proxy(pd$dists$dist_s, pw$Wd)
  
  # Oracle equality
  W_oracle <- oracle_prep_w_GD(pd$dists, H = n)
  expect_equal(pw$Wd, W_oracle, tolerance = 1e-12)
})

test_that("prep_w: GD gauss, kNN (NN<n), oracle equality", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 1:n
  NN <- 15
  H  <- 8
  
  pd <- prep_d(coords = data$coords, NN = NN, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = "gauss", Type = "GD")
  
  pw <- prep_w(H = c(H), kernels = "gauss", Type = "GD",
               adaptive = FALSE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  check_Wd_basic(pw$Wd)
  check_monotone_decreasing_proxy(pd$dists$dist_s, pw$Wd)
  
  W_oracle <- oracle_prep_w_GD(pd$dists, H = H)
  expect_equal(pw$Wd, W_oracle, tolerance = 1e-12)
})

# =============================================================================
# GROUP 2: T (temporal)
# =============================================================================

test_that("prep_w: T gauss, NN=n, oracle equality + monotone proxy", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  
  pd <- prep_d(coords = as.matrix(data$coords_3col[, 3], ncol = 1),
               NN = n, TP = TP, extrapol = FALSE, QP = NULL,
               kernels = "gauss", Type = "T")
  
  pw <- prep_w(H = c(n), kernels = "gauss", Type = "T",
               adaptive = FALSE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  check_Wd_basic(pw$Wd)
  check_monotone_decreasing_proxy(pd$dists$dist_t, pw$Wd)
  
  W_oracle <- oracle_prep_w_T(pd$dists, kernel_name = "gauss", H = n)
  expect_equal(pw$Wd, W_oracle, tolerance = 1e-12)
})

test_that("prep_w: T gauss_SYM_365 uses cyclic distances (oracle equality)", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  
  pd <- prep_d(coords = as.matrix(data$coords_3col[, 3], ncol = 1),
               NN = n, TP = TP, extrapol = FALSE, QP = NULL,
               kernels = "gauss_SYM_365", Type = "T")
  
  pw <- prep_w(H = c(n), kernels = "gauss_SYM_365", Type = "T",
               adaptive = FALSE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  check_Wd_basic(pw$Wd)
  check_monotone_decreasing_proxy(pd$dists$dist_t, pw$Wd)
  
  W_oracle <- oracle_prep_w_T(pd$dists, kernel_name = "gauss", H = n)
  expect_equal(pw$Wd, W_oracle, tolerance = 1e-12)
})

test_that("prep_w: T 'past' does not change weights when dist_t is absolute (same dists)", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  
  # On construit UNE SEULE FOIS les distances, en cyclique par exemple
  pd <- prep_d(
    coords = as.matrix(data$coords_3col[,3], ncol = 1),
    NN = n, TP = TP, extrapol = FALSE, QP = NULL,
    kernels = "gauss", Type = "T"
  )
  
  # Même dists, seul le format change dans prep_w()
  pw_sym  <- prep_w(H = c(n), kernels = "gauss",  Type = "T",
                    adaptive = FALSE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  pw_past <- prep_w(H = c(n), kernels = "gauss_past", Type = "T",
                    adaptive = FALSE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  # Comme dist_t est abs(...) => past=(dist_t>=0) est tout à 1 => identique
  expect_equal(pw_past$Wd, pw_sym$Wd, tolerance = 1e-12)
})
# =============================================================================
# GROUP 3: GDT (spatiotemporal) - alpha mixing
# =============================================================================

test_that("prep_w: GDT alpha=1 equals product oracle; alpha=0 equals sum oracle", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  
  pd <- prep_d(coords = data$coords_3col, NN = n, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = c("gauss","gauss"), Type = "GDT")
  
  pw_prod <- prep_w(H = c(n, n), kernels = c("gauss","gauss"), Type = "GDT",
                    adaptive = c(FALSE,FALSE), dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  pw_sum  <- prep_w(H = c(n, n), kernels = c("gauss","gauss"), Type = "GDT",
                    adaptive = c(FALSE,FALSE), dists = pd$dists, indexG = pd$indexG, alpha = 0)
  
  W_oracle_prod <- oracle_prep_w_GDT(pd$dists, H = c(n,n), alpha = 1)
  W_oracle_sum  <- oracle_prep_w_GDT(pd$dists, H = c(n,n), alpha = 0)
  
  expect_equal(pw_prod$Wd, W_oracle_prod, tolerance = 1e-12)
  expect_equal(pw_sum$Wd,  W_oracle_sum,  tolerance = 1e-12)
  
  check_Wd_basic(pw_prod$Wd)
  check_Wd_basic(pw_sum$Wd)
})

test_that("prep_w: GDT alpha in (0,1) matches mixture oracle", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  a  <- 0.35
  
  pd <- prep_d(coords = data$coords_3col, NN = n, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = c("gauss","gauss"), Type = "GDT")
  
  pw_mix <- prep_w(H = c(n, n), kernels = c("gauss","gauss"), Type = "GDT",
                   adaptive = c(FALSE,FALSE), dists = pd$dists, indexG = pd$indexG, alpha = a)
  
  W_oracle_mix <- oracle_prep_w_GDT(pd$dists, H = c(n,n), alpha = a)
  expect_equal(pw_mix$Wd, W_oracle_mix, tolerance = 1e-12)
  check_Wd_basic(pw_mix$Wd)
})

# =============================================================================
# GROUP 4: GDT cyclic temporal kernel (SYM_365)
# =============================================================================

test_that("prep_w: GDT with cyclic temporal kernel uses cyclic dist_t (oracle equality)", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  
  pd <- prep_d(coords = data$coords_3col, NN = n, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = c("gauss","gauss_SYM_365"), Type = "GDT")
  
  pw <- prep_w(H = c(n, n), kernels = c("gauss","gauss_SYM_365"), Type = "GDT",
               adaptive = c(FALSE,FALSE), dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  # oracle uses gauss on already-cyclic dist_t
  W_oracle <- oracle_prep_w_GDT(pd$dists, H = c(n,n), alpha = 1,
                                kernel_s = "gauss", kernel_t = "gauss")
  expect_equal(pw$Wd, W_oracle, tolerance = 1e-12)
  check_Wd_basic(pw$Wd)
})

# =============================================================================
# GROUP 5: adaptive behavior (basic sanity)
# =============================================================================

test_that("prep_w: GD adaptive=TRUE rounds H and returns valid weights", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 1:n
  NN <- 15
  
  pd <- prep_d(coords = data$coords, NN = NN, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = "gauss", Type = "GD")
  
  pw <- prep_w(H = c(7.7), kernels = "gauss", Type = "GD",
               adaptive = TRUE, dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  check_Wd_basic(pw$Wd)
  expect_equal(dim(pw$Wd), dim(pd$dists$dist_s))
})

test_that("prep_w: GDT adaptive vector accepted and returns valid weights", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  
  pd <- prep_d(coords = data$coords_3col, NN = n, TP = TP,
               extrapol = FALSE, QP = NULL, kernels = c("gauss","gauss"), Type = "GDT")
  
  pw <- prep_w(H = c(10.2, 12.9), kernels = c("gauss","gauss"), Type = "GDT",
               adaptive = c(TRUE, TRUE), dists = pd$dists, indexG = pd$indexG, alpha = 1)
  
  check_Wd_basic(pw$Wd)
  expect_equal(dim(pw$Wd), dim(pd$dists$dist_s))
})

# =============================================================================
# End of test-prep_w.R
# =============================================================================