# Out-of-sample prediction must honour the boundary ("corner") solutions of
# TDS_MGWR when the neighbourhood is capped (NN < n). The largest grid value of
# each axis is then the cap (adaptive kernel) or the largest distance within the
# capped neighbourhoods (fixed kernel), not the whole-sample value, and
# predict(method_pred = 'model') must recognise it:
#   - global component (both sentinels): constant coefficient;
#   - spatial-only component (temporal sentinel): same value at a given place
#     whatever the date;
#   - temporal-only component (spatial sentinel): same value at a given date
#     whatever the place.
# Bandwidths are pinned (control_tds$H / $Ht, Inf = no weighting on that axis)
# so that the test does not depend on the search.

library(testthat)
library(mgwrsar)

corner_data <- function() {
  set.seed(7)
  n <- 400
  u <- runif(n); v <- runif(n); tt <- sample(1:60, n, TRUE)
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  Y <- (1 + 2 * u) + 2 * X1 + sin(tt / 8) * X2 + (u + sin(tt / 10)) * X3 + rnorm(n, sd = 0.5)
  list(dat = data.frame(Y, X1, X2, X3, tt), coords = cbind(u, v))
}

# 4 targets, the same 4 places at other dates, and other places at the same dates
corner_targets <- function() {
  set.seed(8)
  k <- 4
  nd <- data.frame(X1 = rnorm(3 * k), X2 = rnorm(3 * k), X3 = rnorm(3 * k))
  t0 <- c(5, 20, 35, 50)
  list(nd = nd, k = k,
       xy = rbind(cbind(runif(k), runif(k)), matrix(NA, k, 2), cbind(runif(k), runif(k))),
       tt = c(t0, t0 + 7, t0))
}

for (ad in c(TRUE, FALSE)) {
  test_that(sprintf("tds_mgtwr with NN < n predicts corner components consistently (adaptive spatial kernel = %s)", ad), {
    g <- corner_data(); tg <- corner_targets(); k <- tg$k
    xy <- tg$xy; xy[(k + 1):(2 * k), ] <- xy[1:k, ]
    H <- if (ad) c(30, Inf, Inf, 60) else c(0.15, Inf, Inf, 0.25)
    m <- suppressWarnings(TDS_MGWR(
      Y ~ X1 + X2 + X3, Model = "tds_mgwr", Ht = 9999, data = g$dat, coords = g$coords,
      kernels = c("gauss", "gauss"),
      control_tds = list(H = H, Ht = c(Inf, Inf, 8, 15), get_AIC = FALSE, ncore = 1, verbose = FALSE),
      control = list(Z = g$dat$tt, adaptive = c(ad, FALSE), Type = "GDT", NN = 150, ncore = 1,
                     alpha = 1, isgcv = FALSE, verbose = FALSE)))
    # pinned infinite bandwidths are mapped to the capped sentinels
    expect_true(m@V[1] < if (ad) nrow(g$dat) else max(dist(g$coords)))
    expect_equal(unname(m@H[c("X1", "X2")]), rep(m@V[1], 2))
    expect_equal(unname(m@Ht[c("Intercept", "X1")]), rep(m@Vt[1], 2))
    expect_equal(sd(m@Betav[, "X1"]), 0)

    nc <- cbind(xy, tg$tt); colnames(nc) <- c("X", "Y", "time")
    B <- as.matrix(predict(m, newdata = tg$nd, newdata_coords = nc, method_pred = "model",
                           exposant = 8, beta_proj = TRUE)$Beta_proj_out)
    same_place <- 1:k; other_date <- (k + 1):(2 * k); other_place <- (2 * k + 1):(3 * k)
    # global: the estimated constant everywhere
    expect_equal(unname(B[, "X1"]), rep(unname(m@Betav[1, "X1"]), 3 * k), tolerance = 1e-10)
    # spatial-only intercept: independent of the date
    expect_equal(unname(B[same_place, "Intercept"]), unname(B[other_date, "Intercept"]), tolerance = 1e-8)
    # temporal-only X2: independent of the place
    expect_equal(unname(B[same_place, "X2"]), unname(B[other_place, "X2"]), tolerance = 1e-8)
    # the space-time component X3 does depend on both (the checks above are not vacuous)
    expect_gt(max(abs(B[same_place, "X3"] - B[other_date, "X3"])), 1e-4)
    expect_gt(max(abs(B[same_place, "X3"] - B[other_place, "X3"])), 1e-4)
  })

  test_that(sprintf("tds_mgwr with NN < n predicts a global component as a constant (adaptive = %s)", ad), {
    g <- corner_data(); tg <- corner_targets()
    H <- if (ad) c(30, Inf, 50, 80) else c(0.15, Inf, 0.2, 0.3)
    m <- suppressWarnings(TDS_MGWR(
      Y ~ X1 + X2 + X3, Model = "tds_mgwr", data = g$dat, coords = g$coords, kernels = "gauss",
      control_tds = list(H = H, get_AIC = FALSE, ncore = 1, verbose = FALSE),
      control = list(adaptive = ad, Type = "GD", NN = 150, ncore = 1, isgcv = FALSE, verbose = FALSE)))
    expect_equal(unname(m@H["X1"]), m@V[1])
    expect_equal(sd(m@Betav[, "X1"]), 0)
    xy <- tg$xy; xy[is.na(xy)] <- 0.5
    B <- as.matrix(predict(m, newdata = tg$nd, newdata_coords = xy, method_pred = "model",
                           exposant = 8, beta_proj = TRUE)$Beta_proj_out)
    expect_equal(unname(B[, "X1"]), rep(unname(m@Betav[1, "X1"]), nrow(B)), tolerance = 1e-10)
    expect_gt(diff(range(B[, "Intercept"])), 1e-4)
  })
}
