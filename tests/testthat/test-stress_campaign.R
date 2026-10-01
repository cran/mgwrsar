# Adversarial campaign (2026-10): estimator configurations (tools/stress_configs.R)
# fitted on data built to hurt local estimators (tools/stress_generators.R):
# repeated or duplicated locations, points on a line, clusters, isolated point,
# anisotropic or geographic coordinates, tiny samples, constant / collinear /
# locally constant / sparse / badly scaled predictors, heavy tails, outliers,
# heteroskedasticity, constant or perfect response, degenerate time axes.
# Every fit must either succeed with the invariants of stress_invariants()
# (finite coefficients and fit, residuals = Y - fit, leverages in [0, 1],
# tS <= n, finite or saturated AICc, finite SE) or stop with an explicit
# message. A fast subset runs always; the full campaign with RUN_LONG_TESTS.

library(testthat)
library(mgwrsar)
source('../../tools/stress_cases.R')
source('../../tools/stress_configs.R')

stress_gens <- stress_generators()
stress_cfgs <- c(stress_configs(), stress_configs2())
cfg_by_label <- setNames(stress_cfgs, vapply(stress_cfgs, `[[`, "", "label"))

# (generator, config) pairs that must stop with an explicit message
stress_expected_errors <- list(
  constant_predictor = c(mgwr_gauss_fixed = "not identifiable", mgwr_bisq_adapt30_SE = "not identifiable",
                         mgwrsar_0_kc_kv = "not identifiable", sar_2sls = "collinear",
                         tds_mgwr_gauss = "collinear", tds_mgwr_bisq_adapt = "collinear", tds_mgtwr_gauss = "collinear",
                         atds_gwr_gauss = "collinear", tds_mgwr_fixed_vars = "collinear", tds_mgtwr_cyclic = "collinear",
                         predict_tds_mgwr = "collinear", summary_mgwr_fdr = "not identifiable",
                         predict_mgwr_model = "not identifiable", bootstrap_test_B4 = "not identifiable"),
  near_collinear     = c(sar_2sls = "collinear", tds_mgwr_gauss = "collinear", tds_mgwr_bisq_adapt = "collinear",
                         tds_mgtwr_gauss = "collinear", atds_gwr_gauss = "collinear", tds_mgwr_fixed_vars = "collinear",
                         tds_mgtwr_cyclic = "collinear", predict_tds_mgwr = "collinear"),
  exact_collinear    = c(sar_2sls = "collinear", tds_mgwr_gauss = "collinear", tds_mgwr_bisq_adapt = "collinear",
                         tds_mgtwr_gauss = "collinear", atds_gwr_gauss = "collinear", tds_mgwr_fixed_vars = "collinear",
                         tds_mgtwr_cyclic = "collinear", predict_tds_mgwr = "collinear"),
  constant_response  = c(mgwrsar_0_kc_kv = "not identifiable"),
  constant_time      = c(tds_mgtwr_gauss = "temporal variation", tds_mgtwr_cyclic = "temporal variation")
)

stress_check <- function(gn, label) {
  d <- stress_gens[[gn]](); cf <- cfg_by_label[[label]]
  expected <- stress_expected_errors[[gn]][label]
  if (is.null(expected)) expected <- NA_character_
  res <- tryCatch(suppressWarnings(cf$fit(d)), error = function(e) e)
  if (!is.na(expected)) {
    expect_true(inherits(res, "error"), info = paste(gn, label, "should stop"))
    if (inherits(res, "error")) expect_match(conditionMessage(res), expected, info = paste(gn, label))
    return(invisible())
  }
  expect_false(inherits(res, "error"), info = paste(gn, label, if (inherits(res, "error")) conditionMessage(res)))
  if (inherits(res, "error")) return(invisible())
  inv <- if (methods::is(res, "mgwrsar")) stress_invariants(res, d) else if (is.list(res) && !is.null(res$invariants)) res$invariants else character()
  expect_length(inv, 0)
  if (length(inv)) cat("\n", gn, label, ":", paste(inv, collapse = " | "), "\n")
}

fast_gens <- c("repeated_sites", "exact_duplicates", "tiny_n", "constant_predictor", "exact_collinear", "local_constant", "constant_time")
fast_cfgs <- c("gwr_gauss_fixed", "gwr_bisq_adapt30", "gwr_gauss_fixed_SE", "mgwr_gauss_fixed", "gtwr_gauss_fixed", "gtwr_bisq_adapt30", "gtwr_hwb2010", "predict_gwr_model")

for (gn in fast_gens) test_that(paste("stress (fast):", gn), {
  for (label in fast_cfgs) if (isTRUE(cfg_by_label[[label]]$applies(stress_gens[[gn]]()))) stress_check(gn, label)
})

test_that("stress (fast): the space-time bandwidth search ends on duplicated observations", {
  # the golden refinement looped forever here (bracket narrower than 1.3 tolerances)
  t0 <- Sys.time()
  stress_check("exact_duplicates", "search_gdt_gauss_aicc")
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 120)
})

test_that("stress (fast): bad inputs stop with explicit messages", {
  for (cf in Filter(function(c) c$family == "bad", stress_cfgs)) {
    res <- cf$fit(stress_gens$repeated_sites())
    expect_length(res$invariants, 0)
    if (length(res$invariants)) cat("\n", cf$label, ":", res$invariants, "\n")
  }
})

test_that("stress (full campaign)", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "long test")
  for (gn in names(stress_gens)) for (cf in stress_cfgs) {
    if (cf$family == "bad") next
    if (cf$family == "post" && !(gn %in% c("repeated_sites", "exact_duplicates", "tiny_n", "constant_predictor", "near_collinear", "outliers", "lonlat_degrees"))) next
    if (isTRUE(cf$applies(stress_gens[[gn]]()))) stress_check(gn, cf$label)
  }
})
