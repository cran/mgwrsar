# The standard errors of the model returned by TDS_MGWR() must be computed
# from the same state as Betav, Shat, R_k, tS and AICc, whichever sweep is
# returned: the terminal one (default since 1.4.3), the best one by RMSE
# (retained_state = "best", the 1.4.1-1.4.2 behaviour, where that sweep need
# not be the last one), or a polished one.

library(testthat)
library(mgwrsar)

se_dat <- local({
  set.seed(2)  # n = 400: best sweep 8 of 11 with mgwrsar 1.4.2
  n <- 400L
  u <- runif(n); v <- runif(n)
  dat <- data.frame(X2 = rnorm(n), X3 = rnorm(n), X4 = rnorm(n))
  dat$Y <- 3 * (u + v) + 2 * u * sin(4 * v) * dat$X2 +
    4 * u * sin(6 * v) * dat$X3 + 2 * (u - 0.5) * dat$X4 + rnorm(n, sd = 1.5)
  list(dat = dat, co = cbind(u, v), n = n)
})

se_fit <- function(extra = list()) {
  suppressMessages(TDS_MGWR(
    Y ~ X2 + X3 + X4, data = se_dat$dat, coords = se_dat$co, Model = "tds_mgwr",
    kernels = "gauss",
    control_tds = modifyList(list(nns = 10, get_AIC = TRUE, ncore = 1,
                                  init_model = "OLS"), extra),
    control = list(adaptive = TRUE, NN = se_dat$n, Type = "GD", SE = TRUE)))
}

se_check <- function(m) {
  n <- se_dat$n
  sigma2 <- sum(m@residuals^2) / (n - m@tS)
  for (k in colnames(m@sev)) {
    se_k <- sqrt(sigma2 * rowSums(m@R_k[[k]]^2) / m@X[, k]^2)
    expect_equal(unname(m@sev[, k]), unname(se_k), tolerance = 1e-10, info = k)
  }
  expect_equal(m@tS, sum(diag(m@Shat)))
  expect_equal(m@AICc, mgwrsar:::aicc_f(m@residuals, m@tS, n))
  expect_identical(m@Betav, m@HBETA[[length(m@HBETA)]])
}

test_that("tds_mgwr: sev is built from the returned R_k (terminal state)", {
  se_check(se_fit())
})

test_that("tds_mgwr: sev is built from the returned R_k (best state)", {
  se_check(se_fit(list(retained_state = "best")))
})

test_that("tds_mgwr: sev is built from the returned R_k (polished state)", {
  se_check(se_fit(list(polish = 2)))
})
