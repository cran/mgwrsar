# =============================================================================
# Tests for the diagnostic levers of TDS_MGWR() (control_tds$H, $Ht, $V,
# $first_nn, $min_dist_t and control$adaptive[2])
# =============================================================================
# Package: mgwrsar
#
# What is checked:
# - control_tds$H and control_tds$Ht fix the bandwidths of a whole axis: H alone
#   leaves the temporal search free, Ht alone leaves the spatial one free; an
#   axis is either fully fixed or fully searched (no NA, no named subset)
# - `$` partial matching: control_tds$Ht given alone must never be read as
#   control_tds$H (nor min_dist_t as min_dist)
# - the returned slots always describe the returned coefficients: a pinned model
#   is never the OLS start wearing the requested bandwidths
# - a bandwidth at or above the largest one (e.g. Inf) means a global bandwidth
#   instead of silently reusing the model of the previous coefficient
# - the AICc of a pinned model is the one of the returned iteration
# - a user-supplied grid V is order-invariant, and first_nn < n does not crash
# - an adaptive temporal kernel (control$adaptive[2] = TRUE) is searched on a
#   grid of neighbour counts, mirroring the spatial one
# - control_tds$TRUEBETA fills the HRMSE slot: one row per kept iteration
#   (starting model included), the RMSE of each coefficient, their mean, the
#   AICc and the residual RMSE of that iteration
# =============================================================================

library(testthat)
library(mgwrsar)

# Spatio-temporal, non-panel DGP: n distinct locations and dates
lev_n <- 200L
set.seed(1001)
lev_u   <- runif(lev_n); lev_v <- runif(lev_n); lev_day <- runif(lev_n, 1, 365)
lev_b1  <- 3 * (lev_u + lev_v)
lev_b2  <- 4 * sin(2 * pi * lev_day / 365) + 2 * lev_u
lev_b3  <- 4 * lev_u * sin(6 * lev_v)
lev_b4  <- 2 * (lev_u - 0.5)
lev_dat <- data.frame(X2 = rnorm(lev_n), X3 = rnorm(lev_n), X4 = rnorm(lev_n))
lev_dat$Y <- lev_b1 + lev_b2 * lev_dat$X2 + lev_b3 * lev_dat$X3 +
  lev_b4 * lev_dat$X4 + rnorm(lev_n)
lev_co   <- cbind(lev_u, lev_v)
lev_ols  <- sqrt(mean(residuals(lm(Y ~ X2 + X3 + X4, lev_dat))^2))
lev_vars <- c("Intercept", "X2", "X3", "X4")

lev_H  <- c(60, 80, 150, 20)
lev_Ht <- c(50, 30, 190, 120)

lev_fit <- function(extra = list(), Type = "GDT", adaptive = c(TRUE, FALSE),
                    control_extra = list()) {
  ct <- modifyList(list(nns = 8, get_AIC = FALSE, verbose = FALSE, ncore = 1,
                        init_model = "OLS"), extra)
  ctrl <- list(adaptive = adaptive, NN = lev_n, Type = Type)
  if (Type == "GDT") ctrl$Z <- lev_day else ctrl$adaptive <- adaptive[1]
  ctrl <- modifyList(ctrl, control_extra)
  suppressMessages(TDS_MGWR(
    formula = Y ~ X2 + X3 + X4, data = lev_dat, coords = lev_co,
    Model = if (Type == "GDT") "tds_mgtwr" else "tds_mgwr",
    kernels = if (Type == "GDT") c("gauss", "gauss") else "gauss",
    control_tds = ct, control = ctrl))
}

# -----------------------------------------------------------------------------
# 1. Both axes pinned
# -----------------------------------------------------------------------------

test_that("H + Ht: the slots carry the requested bandwidths and describe Betav", {
  m <- lev_fit(list(H = lev_H, Ht = lev_Ht))
  expect_equal(unname(m@H), lev_H)
  expect_equal(unname(m@Ht), lev_Ht)
  expect_identical(names(m@H), lev_vars)
  expect_identical(names(m@Ht), lev_vars)
  # never the OLS start wearing the requested bandwidths
  expect_lt(m@RMSE, lev_ols)
  expect_gt(sd(m@Betav[, "X3"]), 0)
  # the slots are not labels: other bandwidths give other coefficients
  m2 <- lev_fit(list(H = c(60, 80, 40, 20), Ht = lev_Ht))
  expect_gt(max(abs(m@Betav[, "X3"] - m2@Betav[, "X3"])), 1e-6)
})

test_that("a bandwidth at or above the largest one is a global bandwidth", {
  m_ref <- lev_fit(list(H = lev_H, Ht = lev_Ht))
  top_t <- m_ref@max_dist_t
  # 1e6 used to match no branch: the previous coefficient's model was reused,
  # the fit never beat OLS and the OLS start was returned
  m_big <- lev_fit(list(H = lev_H, Ht = c(50, 30, 190, 1e6)))
  m_inf <- lev_fit(list(H = lev_H, Ht = c(50, 30, 190, Inf)))
  m_top <- lev_fit(list(H = lev_H, Ht = c(50, 30, 190, top_t)))
  expect_equal(unname(m_big@Ht), c(50, 30, 190, top_t))
  expect_equal(m_big@Betav, m_top@Betav)
  expect_equal(m_inf@Betav, m_top@Betav)
  expect_lt(m_big@RMSE, lev_ols)
  # same on the spatial axis (adaptive: the largest bandwidth is n)
  m_s <- lev_fit(list(H = c(60, 80, Inf, 20), Ht = lev_Ht))
  expect_equal(unname(m_s@H), c(60, 80, lev_n, 20))
  expect_lt(m_s@RMSE, lev_ols)
})

test_that("the AICc of a pinned model is the one of the returned iteration", {
  # the pinning loop used to overwrite the backfitting iteration counter,
  # shifting the AICc lookup by length(varying) - 1 iterations
  for (Type in c("GDT", "GD")) {
    m <- lev_fit(list(H = lev_H, Ht = if (Type == "GDT") lev_Ht, get_AIC = TRUE),
                 Type = Type)
    expect_equal(m@AICc, mgwrsar:::aicc_f(m@residuals, m@tS, lev_n), info = Type)
    expect_false(any(vapply(m@HBETA[-1], is.null, NA)), info = Type)
  }
})

test_that("standard errors are available with fully pinned bandwidths", {
  m <- lev_fit(list(H = lev_H, Ht = lev_Ht), control_extra = list(SE = TRUE))
  expect_identical(dim(m@sev), c(lev_n, 4L))
  expect_false(anyNA(m@sev))
})

# -----------------------------------------------------------------------------
# 2. One axis left free
# -----------------------------------------------------------------------------

test_that("Ht alone pins the temporal axis and leaves the spatial search free", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  m <- lev_fit(list(Ht = lev_Ht))
  expect_equal(unname(m@Ht), lev_Ht)
  # `$H` partially matches `Ht`: the temporal values used to land in the H slot
  expect_false(isTRUE(all.equal(unname(m@H), lev_Ht)))
  expect_true(all(m@H %in% m@V))
  expect_lt(m@RMSE, lev_ols)
})

test_that("H alone pins the spatial axis and leaves the temporal search free", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  m <- lev_fit(list(H = lev_H))
  expect_equal(unname(m@H), lev_H)
  expect_true(all(m@Ht %in% m@Vt))
  expect_false(all(m@Ht == m@max_dist_t))
  # the former total freeze is spelled H + global temporal bandwidths
  m_frozen <- lev_fit(list(H = lev_H, Ht = rep(Inf, 4)))
  expect_equal(unname(m_frozen@Ht), rep(m_frozen@max_dist_t, 4))
  expect_lte(m@RMSE, m_frozen@RMSE)
})

test_that("H alone also fixes the whole spatial axis of a spatial-only model", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  g <- lev_fit(list(H = lev_H), Type = "GD")
  expect_equal(unname(g@H), lev_H)
  expect_lt(g@RMSE, lev_ols)
})

# -----------------------------------------------------------------------------
# 3. Validation
# -----------------------------------------------------------------------------

test_that("malformed pins are rejected with an explicit message", {
  expect_error(lev_fit(list(H = 100)), "one value per varying coefficient")
  expect_error(lev_fit(list(Ht = c(50, 30))), "one value per varying coefficient")
  # partial pins are not supported: an axis is fully fixed or fully searched
  expect_error(lev_fit(list(H = c(X3 = 150, Intercept = 60))), "unnamed vector")
  expect_error(lev_fit(list(H = c(60, 80, NA, 20))), "must not contain NA")
  expect_error(lev_fit(list(H = c(60, 80, -1, 20))), "must be positive")
  expect_error(lev_fit(list(H = c("a", "b", "c", "d"))), "must be numeric")
  expect_error(lev_fit(list(Ht = lev_Ht), Type = "GD"), "only meaningful when")
})

# -----------------------------------------------------------------------------
# 4. Bandwidth grids
# -----------------------------------------------------------------------------

test_that("a user-supplied grid V is order-invariant", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # an increasing V used to crash the search
  grid <- c(10, 20, 40, 80, 120, 160)
  m_up   <- lev_fit(list(V = grid), Type = "GD")
  m_down <- lev_fit(list(V = rev(grid)), Type = "GD")
  expect_identical(m_up@V, c(lev_n, rev(grid)))
  expect_identical(m_up@V, m_down@V)
  expect_equal(m_up@H, m_down@H)
  expect_equal(m_up@Betav, m_down@Betav)
  expect_error(lev_fit(list(V = c(10, NA, 40)), Type = "GD"), "without NA")
})

test_that("first_nn below n does not crash the adaptive search", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # a coefficient settling on first_nn (the top of V, below max_dist = n) had no
  # candidate above it: "replacement has length zero". Space-time searches hit
  # it systematically, a coefficient being often global in space only.
  m <- lev_fit(list(first_nn = 150L))
  expect_identical(m@V[1], 150)
  expect_true(all(m@H <= lev_n))
  expect_lt(m@RMSE, lev_ols)
  g <- lev_fit(list(first_nn = 150L), Type = "GD")
  expect_true(all(g@H <= lev_n))
})

# Panel DGP: a 6 x 6 lattice observed over T = 5 periods (repeated locations)
pan_N <- 36L; pan_T <- 5L
pan_site <- expand.grid(u = seq(0, 1, length.out = 6), v = seq(0, 1, length.out = 6))
pan_co <- as.matrix(pan_site[rep(seq_len(pan_N), pan_T), ])
pan_t <- rep(seq_len(pan_T), each = pan_N)
set.seed(1002)
pan_dat <- data.frame(X2 = rnorm(pan_N * pan_T), X3 = rnorm(pan_N * pan_T))
pan_dat$Y <- 3 * (pan_co[, 1] + pan_co[, 2]) + (4 * sin(6 * pan_co[, 1]) + pan_t / pan_T) * pan_dat$X2 +
  2 * (pan_co[, 2] - 0.5) * pan_dat$X3 + rnorm(pan_N * pan_T)
pan_d1 <- 0.2   # lattice step = distance to the first site

pan_fit <- function(extra = list(), kernel = "gauss") {
  ct <- modifyList(list(nns = 8, get_AIC = FALSE, verbose = FALSE, ncore = 1,
                        init_model = "OLS"), extra)
  suppressMessages(suppressWarnings(TDS_MGWR(
    formula = Y ~ X2 + X3, data = pan_dat, coords = pan_co, Model = "tds_mgtwr",
    kernels = c(kernel, kernel), control_tds = ct,
    control = list(Z = pan_t, adaptive = c(FALSE, FALSE), NN = pan_N * pan_T, Type = "GDT"))))
}

test_that("in a panel the fixed-kernel spatial grid is extended below the first site", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # the neighbour-count grid stops at the 3rd site (one lattice step): with
  # panel_floor it goes on down to 0.2 * d1 (Gaussian), 0.5 * d1 (bi-square)
  m0 <- pan_fit(list(panel_floor = 1))
  expect_gte(min(m0@V), pan_d1 - 1e-8)
  m <- pan_fit()
  expect_lt(min(m@V), pan_d1)
  expect_gte(min(m@V), 0.2 * pan_d1)
  expect_identical(m@V[seq_along(m0@V)], m0@V)
  expect_true(all(diff(m@V) < 0))
  expect_true(all(m@H >= min(m@V)))
  b <- pan_fit(kernel = "bisq")
  expect_lt(min(b@V), pan_d1)
  expect_gte(min(b@V), 0.5 * pan_d1)
  expect_error(pan_fit(list(panel_floor = 0)), "panel_floor")
  expect_error(pan_fit(list(panel_floor = c(0.2, 0.5))), "panel_floor")
})

test_that("a spatial grid given in distances (V_dist) is taken as is", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  vd <- c(0.05, 0.1, 0.3, 0.6, 1)
  m <- pan_fit(list(V_dist = vd))
  expect_identical(m@V[-1], rev(vd))
  expect_gt(m@V[1], 1)                      # the largest distance is added on top
  expect_true(all(m@H %in% m@V))
  expect_error(pan_fit(list(V_dist = c(0.1, NA))), "V_dist")
})

test_that("min_dist_t alone does not filter the spatial grid", {
  # `$min_dist` partially matches `min_dist_t` when the spatial kernel is fixed
  pins <- list(H = rep(0.3, 4), Ht = lev_Ht)
  m_ref <- lev_fit(pins, adaptive = c(FALSE, FALSE))
  m_mdt <- lev_fit(c(pins, list(min_dist_t = 40)), adaptive = c(FALSE, FALSE))
  expect_identical(m_mdt@V, m_ref@V)
  expect_true(all(m_mdt@Vt > 40))
})

# -----------------------------------------------------------------------------
# 5. Adaptive temporal kernel
# -----------------------------------------------------------------------------

test_that("a temporal-only model accepts an adaptive kernel", {
  # prep_w() suffixed the kernel name twice: "gauss_adapt_sorted_adapt_sorted".
  # TDS needs this model as soon as a coefficient is global in space only.
  m <- MGWRSAR(formula = Y ~ X3, data = lev_dat,
               coords = as.matrix(lev_day, ncol = 1), kernels = "gauss", H = 40,
               Model = "GWR",
               control = list(Type = "T", adaptive = TRUE, NN = lev_n, Z = lev_day))
  expect_gt(sd(m@Betav[, "X3"]), 0)
})

test_that("an adaptive temporal kernel is searched on a grid of neighbour counts", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # used to stop on a NULL @Ht: no temporal grid was built when adaptive[2]
  m <- lev_fit(adaptive = c(TRUE, TRUE))
  expect_identical(m@Vt, m@V)            # mirror of the spatial sequence
  expect_equal(m@max_dist_t, lev_n)
  expect_true(all(m@Ht %in% m@Vt))
  expect_lt(m@RMSE, lev_ols)
  # pins are neighbour counts: rounded, and capped at the global bandwidth
  p <- lev_fit(list(H = lev_H, Ht = c(40.4, 25, 1e6, 120)), adaptive = c(TRUE, TRUE))
  expect_equal(unname(p@Ht), c(40, 25, lev_n, 120))
  expect_lt(p@RMSE, lev_ols)
})

# -----------------------------------------------------------------------------
# 6. Pins combined with the other search options
# -----------------------------------------------------------------------------

test_that("pins hold with a GWR or tds_mgwr starting model", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # init_model = 'GWR' read `$model` from golden_search_bandwidth(), which
  # returns `$best_model`: it failed whatever the levers
  for (init in c("GWR", "tds_mgwr")) {
    m <- suppressWarnings(lev_fit(list(init_model = init, H = lev_H, Ht = lev_Ht)))
    expect_equal(unname(m@H), lev_H, info = init)
    expect_equal(unname(m@Ht), lev_Ht, info = init)
    expect_lt(m@RMSE, lev_ols)
  }
})

test_that("the golden-section refinement leaves pinned coefficients alone", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # the refinement used to flag convergence on the last varying coefficient,
  # which a pinned coefficient never reaches
  g <- suppressWarnings(lev_fit(list(refine = TRUE, check_pairs = TRUE,
                                     H = c(60, 80, 150, 40)), Type = "GD"))
  expect_equal(unname(g@H), c(60, 80, 150, 40))
  expect_lt(g@RMSE, lev_ols)
  # an off-grid pin is not dropped by the check_pairs filter
  m <- lev_fit(list(check_pairs = TRUE, H = c(61, 83, 149, 9)))
  expect_equal(unname(m@H), c(61, 83, 149, 9))
  expect_true(all(m@Ht %in% m@Vt))
})

# -----------------------------------------------------------------------------
# 7. State of the backfitting after a rejected sweep
# -----------------------------------------------------------------------------

lev_X <- cbind(1, as.matrix(lev_dat[, c("X2", "X3", "X4")]))
lev_rmse_path <- function(m)
  vapply(m@HBETA[-1], function(b) sqrt(mean((lev_dat$Y - rowSums(b * lev_X))^2)), 0)

test_that("the descent continues after a rejected sweep and returns the best one", {
  # Restarting from the best state after a rejection repeated the rejected
  # sweep (a sweep is deterministic) and froze the descent; 1.3.2 restarted
  # from an inconsistent state. The descent now continues from the rejected
  # state: later sweeps differ from the rejected one, kept sweeps are the
  # running minima of the RMSE path, and the returned model is the best sweep.
  out <- capture.output(m <- lev_fit(list(Ht = lev_Ht, verbose = TRUE)))
  ln <- grep("rmse = ", out, value = TRUE)
  rm_sweep <- as.numeric(sub(".*rmse = ([0-9.]+).*", "\\1", ln))
  kept <- grepl("\\*$", ln)
  rejected <- which(!kept)
  expect_gte(length(rejected), 1)
  for (r in rejected[rejected < length(rm_sweep)])
    expect_false(isTRUE(all.equal(rm_sweep[r + 1], rm_sweep[r])), info = r)
  expect_equal(rm_sweep[kept], cummin(rm_sweep)[kept])
  expect_equal(m@RMSE, min(rm_sweep), tolerance = 1e-6)
})

test_that("fully pinned bandwidths iterate the backfitting to its fixed point", {
  # nothing is searched: the first sweep (or the last one improving the RMSE)
  # used to be returned, a model that depended on the starting coefficients
  m <- lev_fit(list(H = lev_H, Ht = lev_Ht, tol = 1e-6, maxit = 30))
  path <- lev_rmse_path(m)
  expect_gte(length(path), 5)
  expect_identical(m@Betav, m@HBETA[[length(m@HBETA)]])
  expect_equal(m@RMSE, tail(path, 1))
  expect_lt(abs(diff(tail(path, 2))) / tail(path, 1), 1e-6)
  # the same fixed point from another starting model
  g <- suppressWarnings(lev_fit(list(H = lev_H, Ht = lev_Ht, tol = 1e-6,
                                     maxit = 30, init_model = "GWR")))
  expect_lt(max(abs(g@Betav - m@Betav)), 5e-5)
})

# -----------------------------------------------------------------------------
# 8. True coefficients: the HRMSE history
# -----------------------------------------------------------------------------

lev_TB <- cbind(lev_b1, lev_b2, lev_b3, lev_b4)
lev_beta_rmse <- function(b) unname(sqrt(colMeans((lev_TB - b)^2)))
lev_check_hrmse <- function(m, K = 4L) {
  # one row per kept iteration, starting model included
  expect_identical(nrow(m@HRMSE), length(m@HBETA))
  expect_identical(colnames(m@HRMSE),
                   c(paste0("RMSE_", lev_vars), "meanRMSE", "h", "AICc", "RMSE"))
  expect_false(anyNA(m@HRMSE[, 1:(K + 1)]))
  # column K + 1 is the mean of the K coefficient RMSE
  expect_equal(unname(m@HRMSE[, K + 1]), unname(rowMeans(m@HRMSE[, 1:K])))
  # every row describes the coefficients kept at that iteration, the last one
  # the returned model
  for (i in seq_along(m@HBETA)[-1])
    expect_equal(unname(m@HRMSE[i, 1:K]), lev_beta_rmse(m@HBETA[[i]]), info = i)
  last <- nrow(m@HRMSE)
  expect_equal(unname(m@HRMSE[last, 1:K]), lev_beta_rmse(m@Betav))
  expect_equal(unname(m@HRMSE[last, "RMSE"]), m@RMSE)
}

test_that("TRUEBETA fills the HRMSE history, one row per kept iteration", {
  # HRMSE was written in the backfitting loop but never initialized:
  # "object 'HRMSE' not found"
  m <- lev_fit(list(H = lev_H, Ht = lev_Ht, get_AIC = TRUE, TRUEBETA = lev_TB))
  lev_check_hrmse(m)
  expect_gte(nrow(m@HRMSE), 3)
  expect_equal(unname(m@HRMSE[nrow(m@HRMSE), "AICc"]), m@AICc)
  # row 1 is the starting model (OLS)
  ols <- matrix(coef(lm(Y ~ X2 + X3 + X4, lev_dat)), lev_n, 4, byrow = TRUE)
  expect_equal(unname(m@HRMSE[1, 1:4]), lev_beta_rmse(ols))
  expect_true(is.na(m@HRMSE[1, "h"]))
  # spatial-only model, without AICc, and a data.frame of true coefficients
  g <- lev_fit(list(H = lev_H, TRUEBETA = as.data.frame(lev_TB)), Type = "GD")
  lev_check_hrmse(g)
  expect_true(all(is.na(g@HRMSE[, "AICc"])))
  expect_true(all(g@HRMSE[-1, "h"] %in% c(g@V, g@H)))
})

test_that("a malformed TRUEBETA is rejected with an explicit message", {
  expect_error(lev_fit(list(H = lev_H, Ht = lev_Ht, TRUEBETA = lev_TB[, 1:3])),
               "one column per coefficient")
  expect_error(lev_fit(list(H = lev_H, Ht = lev_Ht, TRUEBETA = lev_TB[-1, ])),
               "one row per observation")
  expect_error(lev_fit(list(H = lev_H, Ht = lev_Ht,
                            TRUEBETA = matrix("a", lev_n, 4))),
               "must be a numeric matrix")
})

test_that("the HRMSE history follows a free search and the atds boosting", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  m <- lev_fit(list(get_AIC = TRUE, TRUEBETA = lev_TB))
  lev_check_hrmse(m)
  expect_equal(unname(m@HRMSE[nrow(m@HRMSE), "AICc"]), m@AICc)
  # stage 2 of atds_mgwr appends its rounds to the stage-1 history
  a <- suppressMessages(TDS_MGWR(
    formula = Y ~ X2 + X3 + X4, data = lev_dat, coords = lev_co,
    Model = "atds_mgwr", kernels = "gauss",
    control_tds = list(nns = 8, get_AIC = TRUE, verbose = FALSE, ncore = 1,
                       init_model = "OLS", TRUEBETA = lev_TB),
    control = list(adaptive = TRUE, NN = lev_n, Type = "GD")))
  expect_gt(nrow(a@HRMSE), 1)
  expect_false(anyNA(a@HRMSE[, 1:5]))
  expect_equal(unname(a@HRMSE[, 5]), unname(rowMeans(a@HRMSE[, 1:4])))
  expect_equal(unname(a@HRMSE[nrow(a@HRMSE), 1:4]), lev_beta_rmse(a@Betav))
})

# -----------------------------------------------------------------------------
# 9. predict() with a temporal bandwidth stored as a single unnamed value
# -----------------------------------------------------------------------------
# When no sweep is kept, or when the starting model comes from a nested
# TDS_MGWR() call, Ht could be a scalar or an unnamed vector; predict() reads
# the temporal bandwidths by coefficient name and failed. The prediction must
# be the same as with the named vector (bug report of 2026-10-01, annex).

test_that("predict() accepts a scalar or unnamed Ht", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  m <- lev_fit(list(H = lev_H, Ht = rep(90, 4)))
  expect_named(m@Ht, lev_vars)
  new <- 1:25
  nc <- cbind(lev_co[new, ], lev_day[new]); colnames(nc) <- c("X", "Y", "time")
  p_named <- predict(m, newdata = lev_dat[new, ], newdata_coords = nc,
                     method_pred = "model", beta_proj = TRUE, exposant = 8)$Y_predicted
  m_scalar <- m; m_scalar@Ht <- 90
  m_unnamed <- m; m_unnamed@Ht <- unname(m@Ht)
  for (mm in list(m_scalar, m_unnamed)) {
    p <- predict(mm, newdata = lev_dat[new, ], newdata_coords = nc,
                 method_pred = "model", beta_proj = TRUE, exposant = 8)$Y_predicted
    expect_equal(p, p_named)
  }
})

# -----------------------------------------------------------------------------
# 10. Non-monotone descent after a rejected sweep (bug report of 2026-10-01)
# -----------------------------------------------------------------------------
# Repeated locations, correlated predictors, gaussian space-time kernels with
# fixed bandwidths: sweep 18 is rejected. Restarting from the best state
# repeated it and returned the state of sweep 17 (RMSE 0.5362); continuing
# from the rejected state finds a smaller spatial bandwidth for the intercept
# and converges to RMSE 0.5288, the solution of 1.3.2.

test_that("the descent escapes a rejected sweep (coauthor's seed 5)", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  source('../../tools/singular_cases.R', local = TRUE)
  d <- singular_data_sites(seed = 5)
  m <- suppressMessages(suppressWarnings(TDS_MGWR(
    formula = Y ~ F1 + F2 + F3, Model = "tds_mgwr", data = d,
    coords = as.matrix(d[, c("x", "y")]), kernels = c("gauss", "gauss"),
    fixed_vars = NULL,
    control_tds = list(nns = 25, verbose = FALSE, tol = 1e-4, init_model = "OLS"),
    control = list(Z = d$time, NN = nrow(d), adaptive = c(FALSE, FALSE),
                   Type = "GDT", ncore = 1))))
  expect_equal(unname(round(m@H)), c(40155, 35947, 410504, 410504))
  expect_lt(sqrt(mean((d$Y - m@fit)^2)), 0.530)
})
