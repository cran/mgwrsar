# Local standard errors (@sev) against a hand-computed reference.
#
# A local fit at point i is beta_i = A X'W y, with A = (X'WX)^-1, so that
#   Var(beta_i) = sigma^2 C C', C = A X'W,
# whatever the scale of the weights. sigma^2 is the local residual variance,
#   sum(w e^2) / (sum(w) - tr(A X'W^2 X)),
# whose denominator is the expectation of sum(w e^2) under a homoskedastic
# local model: it equals n - p for equal weights. Without the fix the engines
# returned sum(w e^2) / (m - p) * A, which with normalized weights is about
# sigma^2 / m and ignores the bandwidth when NN = n.
# The data are deterministic: no RNG.

sev_n <- 150L
sev_i <- seq_len(sev_n)
sev_co <- cbind(X = sin(sev_i * 1.37) + 0.3 * cos(sev_i * 0.41), Y = cos(sev_i * 2.11) - 0.2 * sin(sev_i * 0.73))
sev_d <- data.frame(x1 = sin(sev_i * 0.91) + 1.5, x2 = cos(sev_i * 1.73), xc = sin(sev_i * 2.39))
sev_noise <- sin(sev_i^2 * 0.017)
sev_d$Y <- (1 + sev_co[, 1]) * sev_d$x1 + (2 - sev_co[, 2]) * sev_d$x2 + 0.8 * sev_d$xc + 0.5 * sev_noise
sev_D <- as.matrix(dist(sev_co))

sev_local <- function(X, y, w) {
  X <- as.matrix(X); A <- solve(crossprod(X, X * w))
  C <- A %*% t(X * w); b <- drop(C %*% y); e <- drop(y - X %*% b)
  s2 <- sum(w * e^2) / (sum(w) - sum(diag(A %*% crossprod(X, X * w^2))))
  list(beta = b, se = sqrt(s2 * rowSums(C^2)))
}
sev_fit <- function(formula, H, fixed_vars = NULL, Model = "GWR", kernels = "gauss") {
  MGWRSAR(formula = formula, data = sev_d, coords = sev_co, fixed_vars = fixed_vars, kernels = kernels, H = H,
          Model = Model, control = list(SE = TRUE, adaptive = FALSE, NN = sev_n, isgcv = FALSE, verbose = FALSE))
}
sev_w <- function(H) exp(-0.5 * (sev_D / H)^2)

test_that("univariate GWR: @sev matches the hand-computed variance", {
  H <- 0.5; m <- sev_fit(Y ~ -1 + x1, H); W <- sev_w(H)
  ref <- lapply(sev_i, function(i) sev_local(sev_d["x1"], sev_d$Y, W[i, ]))
  expect_equal(as.numeric(m@Betav), unname(sapply(ref, `[[`, "beta")), tolerance = 1e-8)
  expect_equal(as.numeric(m@sev), unname(sapply(ref, `[[`, "se")), tolerance = 1e-8)
})

test_that("multivariate GWR: @sev matches the hand-computed variance", {
  H <- 0.5; m <- sev_fit(Y ~ x1 + x2, H); W <- sev_w(H); X <- cbind(1, sev_d$x1, sev_d$x2)
  ref <- t(sapply(sev_i, function(i) sev_local(X, sev_d$Y, W[i, ])$se))
  expect_equal(unname(m@sev), unname(ref), tolerance = 1e-8)
})

test_that("mixed GWR: the local residual is net of the fixed part", {
  H <- 0.5; m <- sev_fit(Y ~ x1 + x2 + xc, H, fixed_vars = "xc", Model = "MGWR"); W <- sev_w(H)
  X <- cbind(1, sev_d$x1, sev_d$x2); ystar <- sev_d$Y - sev_d$xc * m@Betac[["xc"]]
  ref <- lapply(sev_i, function(i) sev_local(X, ystar, W[i, ]))
  expect_equal(unname(m@Betav), unname(t(sapply(ref, `[[`, "beta"))), tolerance = 1e-8)
  expect_equal(unname(m@sev), unname(t(sapply(ref, `[[`, "se"))), tolerance = 1e-8)
})

test_that("standard errors grow as the bandwidth shrinks, even with NN = n", {
  se_at <- function(H) mean(sev_fit(Y ~ x1 + x2, H)@sev)
  s <- sapply(c(2, 0.8, 0.3), se_at)
  expect_true(all(diff(s) > 0))
  expect_gt(s[3] / s[1], 1.2)
})

test_that("the wide-bandwidth limit is the OLS standard error", {
  ols <- summary(lm(Y ~ x1 + x2, data = sev_d))$coefficients[, "Std. Error"]
  m <- sev_fit(Y ~ x1 + x2, 1e4)
  expect_equal(unname(colMeans(m@sev)), unname(ols), tolerance = 1e-4)
  m1 <- sev_fit(Y ~ -1 + x1, 1e4)
  expect_equal(mean(m1@sev), summary(lm(Y ~ -1 + x1, data = sev_d))$coefficients[, "Std. Error"], tolerance = 1e-4,
               ignore_attr = TRUE)
})

test_that("a rank-deficient local design keeps a zero on the dropped column", {
  d2 <- sev_d; d2$x3 <- 2 * d2$x1
  m <- suppressWarnings(MGWRSAR(formula = Y ~ -1 + x1 + x3, data = d2, coords = sev_co, fixed_vars = NULL, kernels = "gauss",
                                H = 0.5, Model = "GWR", control = list(SE = TRUE, adaptive = FALSE, NN = sev_n, isgcv = FALSE, verbose = FALSE)))
  expect_true(all(m@sev >= 0)); expect_true(all(is.finite(m@sev)))
  expect_gt(mean(rowSums(m@sev == 0) == 1L), 0.9)
})
