# =============================================================================
# Tests for gtwr_HWB2010() - GTWR of Huang, Wu and Barry (2010)
# =============================================================================
# Package: mgwrsar
#
# What is checked:
# - the (h_ST, tau) <-> (h_S, h_T) change of variables, including its two
#   degenerate branches
# - that the estimator is exactly the Type='GDT' gaussian fixed-bandwidth
#   model at the converted bandwidths (no silent divergence introduced by
#   pre-computing the distances)
# - that tau = 0 collapses onto a plain GWR
# - the locked control keys and the argument validation
# - out-of-sample prediction through the inherited predict method
# =============================================================================

library(testthat)
library(mgwrsar)

data(mydata)

gt_coords <- as.matrix(mydata[, c("x", "y")])
gt_time   <- rep(1:10, length.out = nrow(mydata))
gt_form   <- 'Y_gwr ~ X1 + X2 + X3'
# coords are in metres (diagonal ~ 28 km)
gt_h_st   <- 4000
gt_tau    <- 1e6
# no kNN screening: the pre-selection metric does not depend on tau, so
# screening would make the reference comparisons below inexact (see the
# dedicated test at the end)
gt_NN     <- nrow(mydata)

hwb_to_gdt <- bw_hwb2gdt
gdt_to_hwb <- bw_gdt2hwb

# -----------------------------------------------------------------------------
# 1. Change of variables
# -----------------------------------------------------------------------------

test_that("hwb_to_gdt applies the sqrt(2) kernel convention", {
  # mgwrsar gauss is exp(-.5 (d/h)^2), Huang is exp(-d^2/h_ST^2)
  h <- hwb_to_gdt(1, 1)
  expect_equal(h[1], 1 / sqrt(2))
  expect_equal(h[2], 1 / sqrt(2))
  expect_equal(hwb_to_gdt(2, 4)[2], (2 / sqrt(2)) / 2)
})

test_that("the two parameterisations round-trip exactly", {
  cases <- list(c(1, 1), c(0.12, 0.005), c(5, 100), c(4000, 1e6),
                c(2, 0), c(2, Inf))
  for (p in cases) {
    h <- hwb_to_gdt(p[1], p[2])
    expect_equal(as.numeric(gdt_to_hwb(h[1], h[2])), as.numeric(p),
                 info = paste("h_st =", p[1], "tau =", p[2]))
  }
})

test_that("the degenerate branches follow Huang's equation 10", {
  expect_equal(hwb_to_gdt(2, 0),   c(2 / sqrt(2), Inf))   # mu = 0   -> GWR
  expect_equal(hwb_to_gdt(2, Inf), c(Inf, 2 / sqrt(2)))   # lambda=0 -> TWR
})

# -----------------------------------------------------------------------------
# 2. Equivalence with the Type='GDT' model it reparameterises
# -----------------------------------------------------------------------------

test_that("gtwr_HWB2010 is MGWRSAR(Type='GDT') at the converted bandwidths", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  m <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                    h_st = gt_h_st, tau = gt_tau,
                    control = list(NN = gt_NN, SE = TRUE))

  H <- hwb_to_gdt(gt_h_st, gt_tau)
  ref <- MGWRSAR(gt_form, data = mydata, coords = gt_coords,
                 kernels = c('gauss', 'gauss'), H = H, Model = 'GWR',
                 control = list(Type = 'GDT', Z = gt_time,
                                adaptive = c(FALSE, FALSE),
                                NN = gt_NN, SE = TRUE))

  expect_equal(m@Betav, ref@Betav)
  expect_equal(m@fit,   ref@fit)
  expect_equal(m@AICc,  ref@AICc)
  expect_equal(m@RMSE,  ref@RMSE)

  # bandwidths are stored in both parameterisations and agree
  expect_equal(m@H[1], H[1])
  expect_equal(m@Ht[1], H[2])
  expect_equal(as.numeric(gdt_to_hwb(m@H[1], m@Ht[1])), c(gt_h_st, gt_tau))
})

test_that("coefficients actually vary in space and time", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  m <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                    h_st = gt_h_st, tau = gt_tau, control = list(NN = gt_NN))
  expect_true(all(apply(m@Betav, 2, sd) > 0))
  expect_true(all(is.finite(m@Betav)))
})

test_that("tau = 0 collapses onto a plain GWR", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  m0 <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                     h_st = gt_h_st, tau = 0, control = list(NN = gt_NN))
  g0 <- MGWRSAR(gt_form, data = mydata, coords = gt_coords,
                kernels = 'gauss', H = gt_h_st / sqrt(2), Model = 'GWR',
                control = list(NN = gt_NN))
  expect_equal(m0@Betav, g0@Betav)
})

# -----------------------------------------------------------------------------
# 3. Class, slots and summary
# -----------------------------------------------------------------------------

test_that("the returned object is a gtwr extending mgwrsar", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  m <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                    h_st = gt_h_st, tau = gt_tau, control = list(NN = gt_NN))
  expect_s4_class(m, "gtwr")
  expect_true(is(m, "mgwrsar"))
  expect_equal(m@h_st, gt_h_st)
  expect_equal(m@tau, gt_tau)
  expect_false(m@causal)
  expect_true(is.finite(m@w_tail) && m@w_tail >= 0 && m@w_tail <= 1)
  expect_equal(m@Type, "GDT")
  expect_equal(m@kernels, c("gauss", "gauss"))
  expect_false(any(m@adaptive))

  out <- capture.output(summary(m))
  expect_true(any(grepl("Huang, Wu and Barry", out)))
  expect_true(any(grepl("h_ST", out)))
  expect_true(any(grepl("symmetric", out)))
})

test_that("causal = TRUE switches to the past-only temporal kernel", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  m <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                    h_st = gt_h_st, tau = gt_tau,
                    control = list(NN = gt_NN, causal = TRUE))
  expect_true(m@causal)
  expect_equal(m@kernels, c("gauss", "gauss_past"))
  out <- capture.output(summary(m))
  expect_true(any(grepl("past-only", out)))
})

# -----------------------------------------------------------------------------
# 4. Contract of the control list and of the arguments
# -----------------------------------------------------------------------------

test_that("locked control keys are refused", {
  for (k in c("Type", "kernels", "adaptive", "alpha", "Z")) {
    ctl <- list(NN = 50); ctl[[k]] <- "whatever"
    expect_error(
      gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                   h_st = gt_h_st, tau = gt_tau, control = ctl),
      "sets these control keys itself", info = k)
  }
})

test_that("own control keys do not leak to MGWRSAR", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  # assign_control warns on unknown control names; causal/tail_tol/t.units
  # must be consumed before the call
  expect_no_warning(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                 h_st = gt_h_st, tau = gt_tau,
                 control = list(NN = gt_NN, causal = FALSE, tail_tol = 0.5))
  )
})

test_that("arguments are validated", {
  expect_error(
    gtwr_HWB2010(gt_form, data = mydata,
                 coords = cbind(gt_coords, gt_time), time = gt_time,
                 h_st = gt_h_st, tau = gt_tau),
    "exactly 2 columns")
  expect_error(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                 h_st = -1, tau = gt_tau),
    "h_st must be")
  expect_error(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                 h_st = gt_h_st, tau = -1),
    "tau must be")
  expect_error(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords,
                 time = gt_time[-1], h_st = gt_h_st, tau = gt_tau),
    "same length")
  expect_error(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords,
                 time = as.Date("2020-01-01") + gt_time,
                 h_st = gt_h_st, tau = gt_tau),
    "control\\$t.units")
})

test_that("date time indexes are converted with the declared unit", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  m_num <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords,
                        time = gt_time, h_st = gt_h_st, tau = gt_tau,
                        control = list(NN = gt_NN))
  m_dat <- gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords,
                        time = as.Date("2020-01-01") + gt_time - 1,
                        h_st = gt_h_st, tau = gt_tau,
                        control = list(NN = gt_NN, t.units = "days"))
  # same elapsed time in days as the 1..10 integer index
  expect_equal(m_num@Betav, m_dat@Betav)
})

# -----------------------------------------------------------------------------
# 5. kNN screening guard
# -----------------------------------------------------------------------------

test_that("a truncating kNN pre-selection is reported", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  # huge bandwidth + tight screening: the outer ring still carries weight
  expect_warning(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                 h_st = 1e6, tau = 1e6, control = list(NN = 30)),
    "truncate the gaussian tail")
})

test_that("no warning when all observations are retained", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  n <- nrow(mydata)
  expect_no_warning(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                 h_st = 1e6, tau = 1e6, control = list(NN = n))
  )
})

# -----------------------------------------------------------------------------
# 6. Out-of-sample prediction through the inherited method
# -----------------------------------------------------------------------------

test_that("predict works out-of-sample on (x, y, t)", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  idx <- 1:800
  train <- mydata[idx, ]; test <- mydata[-idx, ]
  tr_t  <- gt_time[idx];  te_t <- gt_time[-idx]

  m <- gtwr_HWB2010(gt_form, data = train,
                    coords = as.matrix(train[, c("x", "y")]), time = tr_t,
                    h_st = gt_h_st, tau = gt_tau, control = list(NN = gt_NN))

  p <- predict(m, newdata = test,
               newdata_coords = as.matrix(cbind(test[, c("x", "y")], te_t)),
               method_pred = 'model')
  expect_length(p, nrow(test))
  expect_true(all(is.finite(p)))
  expect_lt(sqrt(mean((p - test$Y_gwr)^2)), sd(test$Y_gwr))
})

# -----------------------------------------------------------------------------
# 7. The kNN pre-selection does not depend on tau
# -----------------------------------------------------------------------------

test_that("with screening on, tau = 0 no longer matches a plain GWR", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  # prep_d ranks neighbours in a joint space-time metric standardised by
  # sample sd, whatever tau is. With NN < n the retained set therefore differs
  # from the NN spatially nearest neighbours a Type='GD' GWR would use, and the
  # tau = 0 identity only holds without screening. This is documented
  # behaviour, not a bug: it is exactly what w_tail is there to flag.
  NNs <- 300
  m0 <- suppressWarnings(
    gtwr_HWB2010(gt_form, data = mydata, coords = gt_coords, time = gt_time,
                 h_st = gt_h_st, tau = 0, control = list(NN = NNs)))
  g0 <- MGWRSAR(gt_form, data = mydata, coords = gt_coords,
                kernels = 'gauss', H = gt_h_st / sqrt(2), Model = 'GWR',
                control = list(NN = NNs))
  expect_false(isTRUE(all.equal(m0@Betav, g0@Betav)))

  # and the guard fires, because the gaussian has not decayed by the boundary
  expect_gt(m0@w_tail, 0.01)
})

# -----------------------------------------------------------------------------
# 8. Calibration through search_bandwidths, the entry point of the package
# -----------------------------------------------------------------------------
# There is no GTWR-specific bandwidth search: the (h_ST, tau) pair is converted
# to (h_S, h_T), search_bandwidths does the work, and the result is read back.
# The setting is the one of the vignette "GWR and MGWR with Space-Time Kernels":
# normalised coordinates and a time index in [0, 1].

bw_setup <- function(n = 300) {
  s <- simu_multiscale(n = n, myseed = 1, type = 'GG2024', constant = NULL,
                       nuls = NULL, config_beta = 'spatiotemp_old',
                       config_snr = 0.9, config_eps = 'normal')
  list(data = s$mydata, coords = s$coords, time = s$mydata$time,
       formula = as.formula('Y ~ X1 + X2 + X3'))
}

test_that("bw_hwb2gdt turns a (h_st, tau) pair into search_bandwidths ranges", {
  expect_equal(bw_hwb2gdt(0.06 * sqrt(2), (0.06 / 0.4)^2), c(0.06, 0.4))
  expect_equal(bw_gdt2hwb(0.06, 0.4), c(0.06 * sqrt(2), (0.06 / 0.4)^2))
  expect_equal(bw_gdt2hwb(bw_hwb2gdt(0.3, 7)[1], bw_hwb2gdt(0.3, 7)[2]),
               c(0.3, 7))
})

test_that("as_gtwr reads a calibrated GDT model in Huang's parameterisation", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  s <- bw_setup()
  res <- search_bandwidths(
    formula = s$formula, data = s$data, coords = s$coords, fixed_vars = NULL,
    kernels = c('gauss', 'gauss'), Model = 'GWR',
    control = list(Z = s$time, criterion = 'AICc', Type = 'GDT',
                   adaptive = c(FALSE, FALSE), alpha = 1),
    hs_range = c(0, 1.4), ht_range = c(0, 1),
    n_seq = 5, n_rounds = 1, refine = FALSE, ncore = 1)

  m <- as_gtwr(res$best_model)
  expect_s4_class(m, "gtwr")
  expect_true(is(m, "mgwrsar"))
  expect_equal(c(m@h_st, m@tau),
               as.numeric(bw_gdt2hwb(res$minimum[1], res$minimum[2])))
  expect_equal(c(m@H[1], m@Ht[1]), as.numeric(res$minimum))
  expect_false(m@causal)
  expect_true(is.na(m@w_tail))
  # nothing but the parameterisation changed
  expect_equal(m@Betav, res$best_model@Betav)
  expect_equal(m@AICc,  res$best_model@AICc)

  out <- capture.output(summary(m))
  expect_true(any(grepl("Huang, Wu and Barry", out)))
  expect_true(any(grepl("h_ST", out)))
})

test_that("refitting with gtwr_HWB2010 at the calibrated pair is identical", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  s <- bw_setup()
  res <- search_bandwidths(
    formula = s$formula, data = s$data, coords = s$coords, fixed_vars = NULL,
    kernels = c('gauss', 'gauss'), Model = 'GWR',
    control = list(Z = s$time, criterion = 'AICc', Type = 'GDT',
                   adaptive = c(FALSE, FALSE), alpha = 1),
    hs_range = c(0, 1.4), ht_range = c(0, 1),
    n_seq = 5, n_rounds = 1, refine = FALSE, ncore = 1)

  hw <- bw_gdt2hwb(res$minimum[1], res$minimum[2])
  m  <- gtwr_HWB2010(s$formula, data = s$data, coords = s$coords,
                     time = s$time, h_st = hw[1], tau = hw[2])

  expect_equal(m@Betav, res$best_model@Betav)
  expect_equal(m@AICc,  res$best_model@AICc)
  expect_equal(m@RMSE,  res$best_model@RMSE)
  # and it carries the screening diagnostic that a coerced model cannot
  expect_true(is.finite(m@w_tail))
})

test_that("as_gtwr refuses models outside the HWB2010 contract", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  s <- bw_setup()
  gd <- MGWRSAR(s$formula, data = s$data, coords = s$coords,
                kernels = 'gauss', H = 0.1, Model = 'GWR',
                control = list(NN = nrow(s$data)))
  expect_error(as_gtwr(gd), "Type = 'GDT'")

  bi <- MGWRSAR(s$formula, data = s$data, coords = s$coords,
                kernels = c('gauss', 'gauss'), H = c(0.1, 0.4), Model = 'GWR',
                control = list(Type = 'GDT', Z = s$time,
                               adaptive = c(FALSE, FALSE), alpha = 0,
                               NN = nrow(s$data)))
  expect_error(as_gtwr(bi), "alpha = 1")

  expect_error(as_gtwr(42), "class mgwrsar")
})

test_that("as_gtwr is idempotent on a gtwr object", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")

  s <- bw_setup()
  m <- gtwr_HWB2010(s$formula, data = s$data, coords = s$coords, time = s$time,
                    h_st = 0.0849, tau = 0.0225)
  expect_identical(as_gtwr(m), m)
})
