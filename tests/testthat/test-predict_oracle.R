# Oracle tests of the prediction modes: every mode of predict() is compared
# with an explicit reference computed by hand on a small synthetic sample
# (dense kernel weights, weighted least squares), not with a stored hash. A
# hash (tests/testthat/_pred_hashes.csv) freezes a result, it does not validate
# it: the TP hashes of versions 1.3.1 to 1.4.3 recorded a wrong prediction
# (the control of the training call was discarded, see NEWS 1.4.4). A relative
# guard between modes completes the oracles: no mode may be far behind the
# best one, nor close to a global OLS.

library(testthat)
library(mgwrsar)

or_n <- 400L; or_nt <- 60L
set.seed(4041)
or_u <- runif(or_n + or_nt); or_v <- runif(or_n + or_nt); or_t <- runif(or_n + or_nt, 0, 10)
or_dat <- data.frame(X2 = rnorm(or_n + or_nt), X3 = rnorm(or_n + or_nt))
or_b0 <- 2 + 3 * or_u + 0.2 * or_t; or_b2 <- 4 * sin(3 * or_u) * cos(2 * or_v) + 0.3 * sin(or_t); or_b3 <- 2 * (or_v - 0.5)
or_dat$Y <- or_b0 + or_b2 * or_dat$X2 + or_b3 * or_dat$X3 + rnorm(or_n + or_nt, sd = 0.3)
or_co <- cbind(or_u, or_v); or_tr <- seq_len(or_n); or_te <- or_n + seq_len(or_nt)
or_nd <- or_dat[or_te, ]; or_nc <- or_co[or_te, ]; or_y <- or_dat$Y[or_te]
or_Xnew <- cbind(1, or_dat$X2[or_te], or_dat$X3[or_te])
or_Dst <- as.matrix(dist(rbind(or_nc, or_co[or_tr, ])))[seq_len(or_nt), -seq_len(or_nt)]   # test x train
or_Dtt <- abs(outer(or_t[or_te], or_t[or_tr], "-"))                                             # temporal, test x train
or_rmse <- function(p) sqrt(mean((as.numeric(p) - or_y)^2))
or_norm <- function(W) W / rowSums(W)
or_gauss <- function(D, h) exp(-0.5 * (D / h)^2)
or_hadapt <- function(D, H) apply(D, 1, function(d) sort(d)[H])   # adaptive bandwidth: distance to the H-th training point
or_sharpen <- function(W, g) or_norm((W / apply(W, 1, max))^g)      # .sharpen_rows: scale by the row maximum, power, normalise
or_ols <- or_rmse(predict(lm(Y ~ X2 + X3, or_dat[or_tr, ]), or_nd))

test_that("GWR, spatial kernel: 'model', 'shepard' and 'tWtp_model' match their hand-made references", {
  H <- 30
  m <- suppressWarnings(MGWRSAR(Y ~ X2 + X3, data = or_dat[or_tr, ], coords = or_co[or_tr, ], fixed_vars = NULL, kernels = "gauss", H = H,
                                Model = "GWR", control = list(adaptive = TRUE, Type = "GD", SE = FALSE, verbose = FALSE)))
  B <- m@Betav
  # 'model': Gaussian kernel centred on the target, bandwidth = distance to its H-th training point, row-normalised
  W <- or_norm(or_gauss(or_Dst, or_hadapt(or_Dst, H)))
  p <- predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = "model", beta_proj = TRUE)
  expect_equal(unname(p$Beta_proj_out), unname(W %*% B), tolerance = 1e-6)
  expect_equal(as.numeric(p$Y_predicted), unname(rowSums((W %*% B) * or_Xnew)), tolerance = 1e-6)
  # 'shepard': the k_extra + 2 = 10 nearest training points, weights (d - d_9)^2 / (d d_9), row-normalised
  k_extra <- 8
  Wsh <- t(apply(or_Dst, 1, function(d) { o <- order(d)[seq_len(k_extra + 2)]; dh <- d[o][k_extra + 1]
    w <- numeric(length(d)); w[o] <- (d[o] - dh)^2 / (d[o] * dh); w }))
  Wsh <- or_norm(Wsh)
  p <- predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = "shepard", k_extra = k_extra, beta_proj = TRUE)
  expect_equal(unname(p$Beta_proj_out), unname(Wsh %*% B), tolerance = 1e-6)
  # 'tWtp_model': kernel centred on each training point i over the training and test points together (its NN = n
  # nearest, itself included, bandwidth = its H-th sorted distance, itself at distance 0), rows normalised over the
  # neighbours of i, then the columns of the test points are transposed and normalised again
  Dall <- as.matrix(dist(rbind(or_co[or_tr, ], or_nc)))   # train then test
  Wi <- matrix(0, or_n, or_n + or_nt)
  for (i in seq_len(or_n)) { d <- Dall[i, ]; o <- order(d)[seq_len(m@NN)]; Wi[i, o] <- or_gauss(d[o], sort(d)[H]) }
  Wt <- or_norm(t(or_norm(Wi)[, or_n + seq_len(or_nt)]))
  p <- predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = "tWtp_model", beta_proj = TRUE)
  expect_equal(unname(p$Beta_proj_out), unname(Wt %*% B), tolerance = 1e-6)
  # relative guard: every mode within 1.5 times the best one, and well ahead of a global OLS
  r <- sapply(c("TP", "tWtp_model", "model", "shepard"), function(mp) or_rmse(suppressWarnings(predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = mp))))
  expect_lt(max(r) / min(r), 1.5)
  expect_lt(max(r), 0.6 * or_ols)
})

test_that("GWR, space-time kernel (GDT): 'model' = sharpened product kernel, 'tWtp_model' identical", {
  H <- c(30, 2); g <- 8   # default exponent
  m <- suppressWarnings(MGWRSAR(Y ~ X2 + X3, data = or_dat[or_tr, ], coords = or_co[or_tr, ], fixed_vars = NULL, kernels = c("gauss", "gauss"), H = H,
                                Model = "GWR", control = list(Z = or_t[or_tr], adaptive = c(TRUE, FALSE), Type = "GDT", SE = FALSE, verbose = FALSE)))
  B <- m@Betav
  Ks <- or_norm(or_gauss(or_Dst, or_hadapt(or_Dst, H[1]))); Kt <- or_norm(or_gauss(or_Dtt, H[2]))
  W <- or_sharpen(or_norm(Ks * Kt), g)
  ndc <- cbind(or_nc, or_t[or_te])
  p <- predict(m, newdata = or_nd, newdata_coords = ndc, method_pred = "model", beta_proj = TRUE)
  expect_equal(unname(p$Beta_proj_out), unname(W %*% B), tolerance = 1e-6)
  expect_equal(as.numeric(p$Y_predicted), unname(rowSums((W %*% B) * or_Xnew)), tolerance = 1e-6)
  p2 <- predict(m, newdata = or_nd, newdata_coords = ndc, method_pred = "tWtp_model", beta_proj = TRUE)
  expect_equal(unname(p2$Beta_proj_out), unname(p$Beta_proj_out), tolerance = 1e-8)
  # a larger exponent concentrates the weights: nearer to the nearest training point's coefficients
  p24 <- predict(m, newdata = or_nd, newdata_coords = ndc, method_pred = "model", exposant = 24, beta_proj = TRUE)
  expect_equal(unname(p24$Beta_proj_out), unname(or_sharpen(or_norm(Ks * Kt), 24) %*% B), tolerance = 1e-6)
  expect_lt(or_rmse(p$Y_predicted), 0.6 * or_ols)
})

test_that("tds_mgwr and tds_mgtwr: 'model' = per-coefficient sharpened kernels at the selected bandwidths", {
  g <- 8
  # spatial, adaptive
  m <- suppressMessages(suppressWarnings(TDS_MGWR(Y ~ X2 + X3, data = or_dat[or_tr, ], coords = or_co[or_tr, ], Model = "tds_mgwr", kernels = "gauss",
                                                  control_tds = list(nns = 8, get_AIC = FALSE, verbose = FALSE, ncore = 1, init_model = "OLS"),
                                                  control = list(adaptive = TRUE, NN = or_n, Type = "GD"))))
  p <- predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = "model", beta_proj = TRUE)
  for (k in colnames(m@Betav)) {
    W <- if (m@H[[k]] == m@V[1]) matrix(1 / or_n, or_nt, or_n) else or_sharpen(or_norm(or_gauss(or_Dst, or_hadapt(or_Dst, round(m@H[[k]])))), g)
    expect_equal(unname(p$Beta_proj_out[, k]), unname(as.numeric(W %*% m@Betav[, k])), tolerance = 1e-6, label = paste("tds_mgwr", k))
  }
  expect_lt(or_rmse(p$Y_predicted), 0.6 * or_ols)
  # space-time
  mt <- suppressMessages(suppressWarnings(TDS_MGWR(Y ~ X2 + X3, data = or_dat[or_tr, ], coords = or_co[or_tr, ], Model = "tds_mgtwr", Ht = 0, kernels = c("gauss", "gauss"),
                                                   control_tds = list(nns = 8, get_AIC = FALSE, verbose = FALSE, ncore = 1, init_model = "OLS"),
                                                   control = list(Z = or_t[or_tr], adaptive = c(TRUE, FALSE), NN = or_n, Type = "GDT"))))
  ndc <- cbind(or_nc, or_t[or_te])
  pt <- predict(mt, newdata = or_nd, newdata_coords = ndc, method_pred = "model", beta_proj = TRUE)
  Ht <- mt@Ht; if (length(Ht) == 1L || is.null(names(Ht))) { Ht <- rep_len(Ht, length(mt@H)); names(Ht) <- names(mt@H) }
  for (k in colnames(mt@Betav)) {
    hs <- mt@H[[k]]; ht <- Ht[[k]]
    Ks <- if (hs == mt@V[1]) matrix(1, or_nt, or_n) else or_norm(or_gauss(or_Dst, or_hadapt(or_Dst, round(hs))))
    Kt <- if (ht == mt@Vt[1]) matrix(1, or_nt, or_n) else or_norm(or_gauss(or_Dtt, ht))
    W <- if (hs == mt@V[1] && ht == mt@Vt[1]) matrix(1 / or_n, or_nt, or_n) else or_sharpen(or_norm(Ks * Kt), g)
    expect_equal(unname(pt$Beta_proj_out[, k]), unname(as.numeric(W %*% mt@Betav[, k])), tolerance = 1e-6, label = paste("tds_mgtwr", k))
  }
  expect_lt(or_rmse(pt$Y_predicted), 0.6 * or_ols)
})

test_that("MGWR (mixed GWR): 'TP' = training fixed coefficients + local fit of the partial residual at the targets", {
  H <- 30
  m <- suppressWarnings(MGWRSAR(Y ~ X2 + X3, data = or_dat[or_tr, ], coords = or_co[or_tr, ], fixed_vars = "X3", kernels = "gauss", H = H,
                                Model = "MGWR", control = list(adaptive = TRUE, Type = "GD", SE = FALSE, verbose = FALSE)))
  # reference: Y - X3 * betac regressed locally on (1, X2) at each target, Gaussian weights, bandwidth = H-th training point
  r_part <- as.numeric(m@Y) - as.numeric(m@XC[, "X3"] * m@Betac[["X3"]])
  XV <- cbind(1, or_dat$X2[or_tr])
  Bv <- t(sapply(seq_len(or_nt), function(i) { h <- sort(or_Dst[i, ])[H]; w <- or_gauss(or_Dst[i, ], h); lm.wfit(XV, r_part, w)$coefficients }))
  p_ref <- or_dat$X3[or_te] * m@Betac[["X3"]] + rowSums(Bv * cbind(1, or_dat$X2[or_te]))
  p <- suppressWarnings(predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = "TP", beta_proj = TRUE))
  expect_equal(as.numeric(p$Y_predicted), p_ref, tolerance = 1e-6)
  expect_equal(unname(p$Beta_proj_out[, c("Intercept", "X2")]), unname(Bv), tolerance = 1e-6)
  expect_equal(unname(p$Beta_proj_out[, "X3"]), rep(m@Betac[["X3"]], or_nt))
  # relative guard between the three modes, and against a global OLS (the mixed model is misspecified here, X3 having a
  # varying coefficient in the DGP: the guard is looser than for GWR)
  r <- sapply(c("TP", "tWtp_model", "model"), function(mp) or_rmse(suppressWarnings(predict(m, newdata = or_nd, newdata_coords = or_nc, method_pred = mp))))
  expect_lt(max(r) / min(r), 1.5)
  expect_lt(max(r), 0.7 * or_ols)
})
