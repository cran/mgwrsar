# predict(method_pred = 'TP') re-estimates the model at the target points with
# the control of the training call. Versions 1.3.1 to 1.4.3 discarded that
# control (is.list() tested on the unevaluated call) and re-estimated with the
# defaults, so an adaptive bandwidth was read as a distance: on unit-square
# coordinates every weight was about 1 and the prediction was a global OLS.

library(testthat)
library(mgwrsar)

tp_n <- 1200L; tp_nt <- 200L
set.seed(2026)
tp_u <- runif(tp_n + tp_nt); tp_v <- runif(tp_n + tp_nt)
tp_dat <- data.frame(X2 = rnorm(tp_n + tp_nt), X3 = rnorm(tp_n + tp_nt))
tp_dat$Y <- (2 + 3 * tp_u) + 4 * sin(3 * tp_u) * cos(2 * tp_v) * tp_dat$X2 + 2 * (tp_v - 0.5) * tp_dat$X3 + rnorm(tp_n + tp_nt, sd = 0.5)
tp_co <- cbind(tp_u, tp_v); tp_tr <- seq_len(tp_n); tp_te <- tp_n + seq_len(tp_nt)
tp_rmse <- function(p) sqrt(mean((as.numeric(p) - tp_dat$Y[tp_te])^2))

test_that("TP prediction of an adaptive GWR keeps the adaptive bandwidth of the training call", {
  m <- suppressWarnings(MGWRSAR(Y ~ X2 + X3, data = tp_dat[tp_tr, ], coords = tp_co[tp_tr, ], fixed_vars = NULL, kernels = "gauss", H = 40,
                                Model = "GWR", control = list(adaptive = TRUE, Type = "GD", SE = FALSE, verbose = FALSE)))
  p_tp <- suppressWarnings(predict(m, newdata = tp_dat[tp_te, ], newdata_coords = tp_co[tp_te, ], method_pred = "TP"))
  # local fit at each target by hand: Gaussian weights, bandwidth = distance to the 40th training point
  D <- as.matrix(dist(rbind(tp_co[tp_te, ], tp_co[tp_tr, ])))[seq_len(tp_nt), -seq_len(tp_nt)]
  X <- cbind(1, tp_dat$X2[tp_tr], tp_dat$X3[tp_tr])
  p_man <- sapply(seq_len(tp_nt), function(i) { h <- sort(D[i, ])[40]; w <- exp(-0.5 * (D[i, ] / h)^2)
    sum(lm.wfit(X, tp_dat$Y[tp_tr], w)$coefficients * c(1, tp_dat$X2[tp_te[i]], tp_dat$X3[tp_te[i]])) })
  expect_equal(as.numeric(p_tp), p_man, tolerance = 1e-6)
  # far better than the global OLS the defaults used to give
  p_ols <- predict(lm(Y ~ X2 + X3, tp_dat[tp_tr, ]), tp_dat[tp_te, ])
  expect_lt(tp_rmse(p_tp), 0.5 * tp_rmse(p_ols))
  # a control passed through a variable (a symbol in the call) is handled too
  ctl <- list(adaptive = TRUE, Type = "GD", SE = FALSE, verbose = FALSE)
  m2 <- suppressWarnings(MGWRSAR(Y ~ X2 + X3, data = tp_dat[tp_tr, ], coords = tp_co[tp_tr, ], fixed_vars = NULL, kernels = "gauss", H = 40,
                                 Model = "GWR", control = ctl))
  rm(ctl)
  p_tp2 <- suppressWarnings(predict(m2, newdata = tp_dat[tp_te, ], newdata_coords = tp_co[tp_te, ], method_pred = "TP"))
  expect_equal(as.numeric(p_tp2), p_man, tolerance = 1e-6)
})

test_that("TP prediction of a fixed-bandwidth GWR matches the hand-made local fit", {
  m <- suppressWarnings(MGWRSAR(Y ~ X2 + X3, data = tp_dat[tp_tr, ], coords = tp_co[tp_tr, ], fixed_vars = NULL, kernels = "gauss", H = 0.08,
                                Model = "GWR", control = list(adaptive = FALSE, Type = "GD", SE = FALSE, verbose = FALSE)))
  p_tp <- suppressWarnings(predict(m, newdata = tp_dat[tp_te, ], newdata_coords = tp_co[tp_te, ], method_pred = "TP"))
  D <- as.matrix(dist(rbind(tp_co[tp_te, ], tp_co[tp_tr, ])))[seq_len(tp_nt), -seq_len(tp_nt)]
  X <- cbind(1, tp_dat$X2[tp_tr], tp_dat$X3[tp_tr])
  p_man <- sapply(seq_len(tp_nt), function(i) { w <- exp(-0.5 * (D[i, ] / 0.08)^2)
    sum(lm.wfit(X, tp_dat$Y[tp_tr], w)$coefficients * c(1, tp_dat$X2[tp_te[i]], tp_dat$X3[tp_te[i]])) })
  expect_equal(as.numeric(p_tp), p_man, tolerance = 1e-6)
})
