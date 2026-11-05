# =============================================================================
# Tests for control_tds$check_pairs (identifiability floor per covariate)
# =============================================================================
# Package: mgwrsar
#
# What is checked:
# - the floor is the local-variance criterion computed through the neighbour
#   index: prep_w() weights are n x NN in neighbour order, so the local moments
#   cannot be taken as `W %*% X` (wrong columns when NN = n, non-conformable
#   when NN < n)
# - a covariate that never reaches the threshold on the grid (here a 2 x 2
#   corner block under a compact kernel: zero weight from the opposite
#   quadrant) is kept at the global bandwidth with a warning naming it, for
#   Type = 'GD' (slot HKmin) and Type = 'GDT' (slot HKMIN, single sentinel
#   pair) alike; the GDT branch used to leave HKMIN[[k]] NULL and crash in
#   update_opt_st()
# - the intercept takes the smallest floor of the covariates (per spatial
#   bandwidth for 'GDT')
# - everything holds with a neighbour cap control$NN < n
# =============================================================================

library(testthat)
library(mgwrsar)

cp_grid <- expand.grid(x = 1:12, y = 1:12)
cp_n <- nrow(cp_grid)
set.seed(7)
cp_dat <- data.frame(
  X1 = rnorm(cp_n),
  Xs = sqrt((cp_grid$x - 1)^2 + (cp_grid$y - 1)^2) / 12,   # smooth distance field
  Xb = as.numeric(cp_grid$x <= 2 & cp_grid$y <= 2))          # locally constant block
cp_dat$Y <- 1 + 2 * cp_dat$X1 + cp_dat$Xs + cp_dat$Xb + rnorm(cp_n)
cp_co <- as.matrix(cp_grid)
cp_tm <- runif(cp_n, 1, 100)

cp_fit <- function(Type, formula = Y ~ X1 + Xs + Xb, NN = cp_n, extra = list()) {
  ct <- modifyList(list(nns = 6, get_AIC = FALSE, verbose = FALSE, ncore = 1,
                        init_model = "OLS", check_pairs = TRUE), extra)
  ctrl <- list(adaptive = c(TRUE, FALSE), NN = NN, Type = Type)
  if (Type == "GDT") ctrl$Z <- cp_tm else ctrl$adaptive <- TRUE
  suppressMessages(TDS_MGWR(
    formula = formula, data = cp_dat, coords = cp_co,
    Model = if (Type == "GDT") "tds_mgtwr" else "tds_mgwr",
    kernels = if (Type == "GDT") c("bisq", "gauss") else "bisq",
    control_tds = ct, control = ctrl))
}

# collect the warnings of an expression, return its value
cp_warn <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  list(value = v, warnings = w)
}

# Independent re-implementation of the 'GD' floor: first bandwidth of the grid
# at which the 10% quantile over sites of local/global variance exceeds 0.05,
# the local moments gathered site by site through indexG.
cp_floor_gd <- function(x, V, NN) {
  G <- mgwrsar:::prep_d(coords = cp_co, NN = NN, TP = 1:cp_n, kernels = "bisq", Type = "GD")
  vg <- var(x)
  for (v in sort(V)) {
    W <- mgwrsar:::prep_w(H = v, kernels = "bisq", Type = "GD", adaptive = TRUE,
                          dists = G$dists, indexG = G$indexG)$Wd
    ratio <- vapply(seq_len(cp_n), function(i) {
      xi <- x[G$indexG[i, seq_len(ncol(W))]]
      (sum(W[i, ] * xi^2) - sum(W[i, ] * xi)^2) / vg
    }, numeric(1))
    if (quantile(ratio, 0.10, na.rm = TRUE) > 0.05) return(v)
  }
  NA_real_
}

# Same for 'GDT': brute force over every (hs, ht) pair of the grid, the first
# admissible ht per hs (NULL when no pair is admissible).
cp_floor_gdt <- function(x, V, Vt, NN) {
  G <- mgwrsar:::prep_d(coords = cbind(cp_co, cp_tm), NN = NN, TP = 1:cp_n,
                        kernels = c("bisq", "gauss"), Type = "GDT")
  vg <- var(x); out <- NULL
  for (hs in sort(V)) for (ht in sort(Vt)) {
    W <- mgwrsar:::prep_w(H = c(hs, ht), kernels = c("bisq", "gauss"), Type = "GDT",
                          adaptive = c(TRUE, FALSE), dists = G$dists, indexG = G$indexG)$Wd
    ratio <- vapply(seq_len(cp_n), function(i) {
      xi <- x[G$indexG[i, seq_len(ncol(W))]]
      (sum(W[i, ] * xi^2) - sum(W[i, ] * xi)^2) / vg
    }, numeric(1))
    if (quantile(ratio, 0.10, na.rm = TRUE) > 0.05) { out <- rbind(out, c(hs = hs, ht = ht)); break }
  }
  out
}

test_that("GD: a locally constant covariate is kept global with a warning, the others get a floor", {
  r <- cp_warn(cp_fit("GD"))
  m <- r$value
  expect_length(r$warnings, 1)
  expect_match(r$warnings, "check_pairs.*Xb.*kept at the global bandwidth")
  expect_equal(unname(m@HKmin["Xb"]), m@max_dist)
  expect_equal(unname(m@H["Xb"]), m@max_dist)
  # the smooth field has an intermediate floor, honoured by the search
  expect_gt(m@HKmin[["Xs"]], min(m@V))
  expect_lt(m@HKmin[["Xs"]], m@max_dist)
  expect_gte(m@H[["Xs"]], m@HKmin[["Xs"]])
  expect_equal(unname(m@HKmin["Intercept"]), min(m@HKmin[-1]))
  # the floor is the criterion computed through the neighbour index
  expect_equal(unname(m@HKmin["Xs"]), cp_floor_gd(cp_dat$Xs, m@V, cp_n))
  expect_equal(unname(m@HKmin["X1"]), cp_floor_gd(cp_dat$X1, m@V, cp_n))
  # no warning without the degenerate covariate
  expect_length(cp_warn(cp_fit("GD", Y ~ X1 + Xs))$warnings, 0)
})

test_that("GDT: a covariate without admissible pair gets the sentinel pair instead of NULL", {
  r <- cp_warn(cp_fit("GDT"))
  m <- r$value
  expect_length(r$warnings, 1)
  expect_match(r$warnings, "check_pairs.*Xb.*kept at the global bandwidth")
  sentinel <- cbind(hs = m@max_dist, ht = m@max_dist_t)
  expect_equal(m@HKMIN[["Xb"]], sentinel, ignore_attr = FALSE)
  expect_equal(unname(m@H["Xb"]), m@max_dist)
  expect_equal(unname(m@Ht["Xb"]), m@max_dist_t)
  # the smooth field has admissible pairs below the global ones and the
  # selected pair is one of them
  fl <- m@HKMIN[["Xs"]]
  expect_true(any(fl[, "hs"] < m@max_dist))
  row <- match(m@H[["Xs"]], fl[, "hs"])
  expect_false(is.na(row))
  expect_gte(m@Ht[["Xs"]], fl[row, "ht"])
  # the boundaries are those of a brute-force scan of the whole (hs, ht) grid
  expect_equal(unname(fl), unname(cp_floor_gdt(cp_dat$Xs, m@V, m@Vt, cp_n)))
  expect_equal(unname(m@HKMIN[["X1"]]), unname(cp_floor_gdt(cp_dat$X1, m@V, m@Vt, cp_n)))
  expect_null(cp_floor_gdt(cp_dat$Xb, m@V, m@Vt, cp_n))
  # intercept: per spatial bandwidth, the smallest temporal floor of the covariates
  rows <- rbind(m@HKMIN[["X1"]], fl)
  expected <- t(sapply(sort(unique(rows[, "hs"])),
                       function(h) c(hs = h, ht = min(rows[rows[, "hs"] == h, "ht"]))))
  expect_equal(unname(m@HKMIN[["Intercept"]]), unname(expected))
})

test_that("check_pairs works with a neighbour cap NN < n (GD and GDT)", {
  NN <- 60L
  r <- cp_warn(cp_fit("GD", NN = NN))
  m <- r$value
  expect_match(r$warnings, "Xb")
  expect_equal(unname(m@HKmin["Xb"]), m@max_dist)
  expect_equal(unname(m@HKmin["Xs"]), cp_floor_gd(cp_dat$Xs, m@V, NN))
  expect_gte(m@H[["Xs"]], m@HKmin[["Xs"]])
  r <- cp_warn(cp_fit("GDT", NN = NN))
  m <- r$value
  expect_match(r$warnings, "Xb")
  expect_equal(m@HKMIN[["Xb"]], cbind(hs = m@max_dist, ht = m@max_dist_t))
  expect_equal(unname(m@H["Xb"]), m@max_dist)
  expect_true(any(m@HKMIN[["Xs"]][, "hs"] < m@max_dist))
})
