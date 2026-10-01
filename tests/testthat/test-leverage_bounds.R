# Leverages of the local weighted least-squares fits are hat-matrix diagonal
# elements, hence in [0, 1], and their sum tS cannot exceed n. The Q-free path
# of the multivariate engine (TS requested without Shat) reconstructed the
# focal row of Q through two triangular solves on R, which explodes when the
# local system is nearly singular: with a space-time kernel and tiny
# bandwidths on repeated locations, 16 leverages out of 400 exceeded 1 (up to
# 10), tS exceeded n - 1 and the AICc penalty changed sign, so the bandwidth
# search settled on the smallest bandwidths (bug report of 2026-10-01).

source('../../tools/singular_cases.R')
leverage_data <- singular_data_sites

test_that("leverages stay in [0, 1] with a space-time kernel at tiny bandwidths", {
  d <- leverage_data()
  coords <- as.matrix(d[, c("x", "y")])
  n <- nrow(d)
  fit <- function(H)
    MGWRSAR(formula = Y ~ F1 + F2 + F3, data = d, coords = coords,
            fixed_vars = NULL, kernels = c("gauss", "gauss"), H = H,
            Model = "GWR",
            control = list(Type = "GDT", adaptive = c(FALSE, FALSE),
                           Z = d$time, criterion = "AICc"))
  m_small <- fit(c(3000, 95))     # nearly interpolating: tS close to n - 1
  m_large <- fit(c(39000, 745))   # the bandwidths selected by AICc in 1.3.2
  for (m in list(m_small, m_large)) {
    expect_true(all(m@TS >= -1e-8 & m@TS <= 1 + 1e-8))
    expect_lte(m@tS, n)
  }
  expect_true(is.finite(m_large@AICc))
  # the saturated fit is rejected (its AICc may be +Inf), not rewarded
  expect_gt(m_small@AICc, m_large@AICc)
})

test_that("aicc_f is +Inf once the trace reaches n - 1", {
  e <- rnorm(50)
  expect_true(is.finite(mgwrsar:::aicc_f(e, ts = 10, n = 50)))
  expect_identical(mgwrsar:::aicc_f(e, ts = 49, n = 50), Inf)
  expect_identical(mgwrsar:::aicc_f(e, ts = 60, n = 50), Inf)
  expect_identical(mgwrsar:::aicc_f(e, ts = 30, n = 50, pena = 2), Inf)
})
