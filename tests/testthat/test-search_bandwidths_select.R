# Best row of a search_bandwidths() grid round (1.4.2): the minimum of the
# finite scores, without the former "median - 10 MAD" screen.

test_that("a genuine minimum below a flat plateau is kept", {
  # round 1 of a GWR AICc search (n = 5000, fixed Gaussian kernel, grid of 10
  # on [0, max_dist]); the zero bandwidth yields no score. The former screen
  # rejected h = 0.15443 (threshold 27761.9) and returned h = 0.30885.
  res <- data.frame(h_s = c(0, 0.15443, 0.30885, 0.46328, 0.6177, 0.77213,
                            0.92655, 1.08098, 1.2354, 1.38983),
                    h_t = NA_real_,
                    score = c(NA, 27726.10, 28340.63, 28645.53, 28779.97, 28845.82,
                              28882.28, 28904.47, 28918.93, 28966.49),
                    isolated = 0)
  best <- mgwrsar:::.sb_select_best(res)
  expect_equal(best$h_s, 0.15443)
  expect_equal(best$score, 27726.10)
})

test_that("non-finite scores are ignored, and the first row is returned when none is finite", {
  res <- data.frame(h_s = 1:4, h_t = NA_real_, score = c(Inf, NA, 10, 12), isolated = 0)
  expect_equal(mgwrsar:::.sb_select_best(res)$h_s, 3)
  res$score <- c(Inf, NA, -Inf, NaN)
  expect_equal(mgwrsar:::.sb_select_best(res)$h_s, 1)
})
