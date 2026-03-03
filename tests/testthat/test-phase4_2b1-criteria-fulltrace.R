test_that("Phase 4.2.B1: full-trace criteria helpers are available", {
  expect_true(exists(".compute_trace_full", mode = "function"))
  expect_true(exists(".compute_enp_from_trace", mode = "function"))
  expect_true(exists(".compute_AICc_from_rss_enp", mode = "function"))

  a <- .compute_AICc_from_rss_enp(n = 100, rss = 120, enp = 5)
  b <- .compute_AICc_from_rss_enp(n = 100, rss = 120, enp = 10)
  expect_true(is.finite(a))
  expect_true(is.finite(b))
  expect_true(b > a)
})
