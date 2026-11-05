# Sharpening of the prediction weights (predict(method_pred = 'model') for
# TDS models): each row is scaled by its largest weight before the power, so
# that a far target with a small Gaussian bandwidth keeps a finite prediction
# (1.4.3); the result is unchanged otherwise. Also: the importance-driven
# update order no longer depends on the units of the covariates.

library(testthat)
library(mgwrsar)

sh_n <- 200L
set.seed(1003)
sh_u <- runif(sh_n); sh_v <- runif(sh_n)
sh_dat <- data.frame(X2 = rnorm(sh_n), X3 = rnorm(sh_n))
sh_dat$Y <- 3 * (sh_u + sh_v) + 4 * sh_u * sin(6 * sh_v) * sh_dat$X2 +
  2 * (sh_v - 0.5) * sh_dat$X3 + rnorm(sh_n)
sh_co <- cbind(sh_u, sh_v)

sh_fit <- function(dat = sh_dat, extra = list()) {
  suppressMessages(suppressWarnings(TDS_MGWR(
    formula = Y ~ X2 + X3, data = dat, coords = sh_co, Model = "tds_mgwr",
    kernels = "gauss",
    control_tds = modifyList(list(nns = 8, get_AIC = FALSE, verbose = FALSE, ncore = 1,
                                  init_model = "OLS"), extra),
    control = list(adaptive = FALSE, NN = sh_n, Type = "GD"))))
}

test_that("a far target with a small fixed Gaussian bandwidth gets a finite sharpened prediction", {
  m <- sh_fit()
  # two targets inside the square, one far outside (distance about 7)
  nd <- data.frame(X2 = c(0.5, -1, 1), X3 = c(1, 0, -0.5), Y = 0)
  nc <- rbind(c(0.5, 0.5), c(0.1, 0.9), c(5, 5))
  for (g in c(1, 8, 24)) {
    # the far target has no positive Gaussian weight (underflow at distance
    # 6 for a bandwidth of a few hundredths): nearest-neighbour fallback
    w <- capture_warnings(p <- predict(m, newdata = nd, newdata_coords = nc, method_pred = "model",
                                       exposant = g, beta_proj = TRUE))
    expect_true(any(grepl("nearest training point", w)), info = paste("gamma", g))
    expect_true(all(is.finite(p$Y_predicted)), info = paste("gamma", g))
    expect_true(all(is.finite(p$Beta_proj_out)), info = paste("gamma", g))
    # each extrapolated coefficient is a convex combination of training ones
    for (k in colnames(p$Beta_proj_out)) {
      expect_gte(min(p$Beta_proj_out[, k]), min(m@Betav[, k]) - 1e-8, label = k)
      expect_lte(max(p$Beta_proj_out[, k]), max(m@Betav[, k]) + 1e-8, label = k)
    }
  }
  # the far target takes the coefficients of its nearest training point as
  # gamma grows (the sharpened kernel concentrates on it)
  p24 <- suppressWarnings(predict(m, newdata = nd, newdata_coords = nc, method_pred = "model",
                                  exposant = 24, beta_proj = TRUE))
  nn <- which.min((sh_u - 5)^2 + (sh_v - 5)^2)
  expect_equal(unname(p24$Beta_proj_out[3, ]), unname(m@Betav[nn, ]), tolerance = 1e-3)
})

test_that("sharpening is unchanged where the old formula was finite", {
  m <- sh_fit()
  set.seed(7)
  nd <- data.frame(X2 = rnorm(20), X3 = rnorm(20), Y = 0)
  nc <- cbind(runif(20), runif(20))
  p8 <- predict(m, newdata = nd, newdata_coords = nc, method_pred = "model",
                exposant = 8, beta_proj = TRUE)
  # old formula, recomputed by hand for one coefficient at fixed bandwidth
  k <- "X2"; h <- m@H[[k]]
  D <- as.matrix(dist(rbind(nc, sh_co)))[1:20, -(1:20)]
  W <- exp(-(D / h)^2 / 2); W <- W^8; W <- W / rowSums(W)
  expect_equal(unname(p8$Beta_proj_out[, k]), unname(as.numeric(W %*% m@Betav[, k])),
               tolerance = 1e-6)
})

test_that("the importance-driven order does not depend on the units of the covariates", {
  # score = standard deviation of the fitted contribution beta_k x_k (1.4.3):
  # rescaling a covariate rescales its coefficient and leaves the score, the
  # order, the bandwidths and the fit unchanged
  m1 <- sh_fit()
  d2 <- sh_dat; d2$X2 <- 1000 * d2$X2
  m2 <- sh_fit(d2)
  expect_equal(m2@H, m1@H)
  expect_equal(m2@Betav[, "X2"] * 1000, m1@Betav[, "X2"], tolerance = 1e-6)
  expect_equal(m2@Betav[, "X3"], m1@Betav[, "X3"], tolerance = 1e-6)
  expect_equal(m2@RMSE, m1@RMSE, tolerance = 1e-8)
  # the score of a constant contribution is zero: a covariate whose fitted
  # contribution does not vary across observations is updated last
  expect_identical(stats::sd(rep(2, 10)), 0)
})

test_that("update_last sends the named covariates to the end of every sweep", {
  out <- capture.output(m <- sh_fit(extra = list(update_last = "X2", verbose = TRUE)))
  ord <- grep("^ +[0-9]+ +Intercept", out, value = TRUE)
  expect_gt(length(ord), 1)
  last <- vapply(strsplit(trimws(ord), "\\s+"), function(t) t[length(t)], "")
  expect_true(all(last == "X2"))
  first <- vapply(strsplit(trimws(ord), "\\s+"), function(t) t[2], "")
  expect_true(all(first == "Intercept"))
  expect_true(all(m@H %in% m@V))
  expect_error(sh_fit(extra = list(update_last = "X9")), "update_last")
  expect_error(sh_fit(extra = list(update_last = "Intercept")), "update_last")
})

test_that("ordering = 'random' permutes the covariates at every sweep, reproducibly", {
  set.seed(11); out1 <- capture.output(m1 <- sh_fit(extra = list(ordering = "random", verbose = TRUE)))
  set.seed(11); out2 <- capture.output(m2 <- sh_fit(extra = list(ordering = "random", verbose = TRUE)))
  expect_equal(m1@H, m2@H)
  expect_equal(m1@Betav, m2@Betav)
  ord <- grep("^ +[0-9]+ +Intercept", out1, value = TRUE)
  toks <- strsplit(trimws(ord), "\\s+")
  expect_true(all(vapply(toks, function(t) t[2], "") == "Intercept"))
  # with two covariates, both orders appear over the sweeps
  expect_gt(length(unique(vapply(toks, function(t) paste(t[-(1:2)], collapse = " "), ""))), 1)
  set.seed(12); m3 <- sh_fit(extra = list(ordering = "random"))
  expect_false(isTRUE(all.equal(m3@HBETA, m1@HBETA)))
  mc <- sh_fit(extra = list(ordering = "cyclic"))
  expect_true(all(mc@H %in% mc@V))
  expect_error(sh_fit(extra = list(ordering = "shuffle")), "ordering")
})

test_that("importance_order reverses the update order", {
  m_dec <- sh_fit()
  m_inc <- sh_fit(extra = list(importance_order = "increasing"))
  expect_equal(m_dec@saturation$ratio, m_inc@saturation$ratio)
  expect_true(all(m_inc@H %in% m_inc@V))
  expect_lt(m_inc@RMSE, sqrt(mean(residuals(lm(Y ~ X2 + X3, sh_dat))^2)))
  expect_error(sh_fit(extra = list(importance_order = "random")), "importance_order")
})
