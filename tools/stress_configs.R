# Estimator configurations of the stress campaign (see tools/stress_cases.R).
# Bandwidths are set from the data (quantiles of the distances) so that every
# generator is fitted at comparable scales; `applies` excludes meaningless
# pairs (a space-time model on a constant time axis is still fitted: it must
# fail cleanly or work).

stress_coords <- function(d) as.matrix(d[, c("x", "y")])
stress_hdist <- function(d, q = 0.3) { D <- dist(stress_coords(d)); as.numeric(quantile(D[D > 0], q)) }
stress_htime <- function(d, q = 0.3) { D <- dist(d$time); v <- D[D > 0]; if (!length(v)) 1 else as.numeric(quantile(v, q)) }
stress_formula <- Y ~ X1 + X2 + X3

stress_mgwrsar <- function(d, Model = "GWR", kernels = "gauss", H, Type = "GD", adaptive = FALSE,
                           fixed_vars = NULL, SE = FALSE, get_s = FALSE, extra = list()) {
  ctl <- c(list(Type = Type, adaptive = adaptive, SE = SE, get_s = get_s, criterion = "AICc"), extra)
  if (Type %in% c("GDT", "T")) ctl$Z <- d$time
  MGWRSAR(formula = stress_formula, data = d, coords = stress_coords(d), fixed_vars = fixed_vars,
          kernels = kernels, H = H, Model = Model, control = ctl)
}

stress_predict <- function(d, fit_fun, frac = 0.25, ...) {
  n <- nrow(d); test <- with_mt_seed(77, sort(sample(n, round(frac * n)))); app <- setdiff(seq_len(n), test)
  m <- fit_fun(d[app, ])
  nc <- cbind(stress_coords(d)[test, , drop = FALSE], d$time[test]); colnames(nc) <- c("X", "Y", "time")
  p <- predict(m, newdata = d[test, ], newdata_coords = nc, ...)
  y <- if (is.list(p) && !is.null(p$Y_predicted)) p$Y_predicted else p
  inv <- character()
  if (length(y) != length(test)) inv <- c(inv, "prediction length != n_test")
  if (!all(is.finite(y))) inv <- c(inv, "prediction not finite")
  if (length(y) == length(test) && all(is.finite(y))) {
    rmse_p <- sqrt(mean((y - d$Y[test])^2)); rmse_lm <- sqrt(mean((predict(lm(stress_formula, d[app, ]), d[test, ]) - d$Y[test])^2))
    if (rmse_p > 10 * rmse_lm + 1e-8) inv <- c(inv, sprintf("prediction RMSE %.3g > 10 x lm %.3g", rmse_p, rmse_lm))
  }
  list(invariants = inv, model = m)
}

stress_configs <- function() list(
  # --- MGWRSAR with fixed bandwidths ---
  stress_config("gwr_gauss_fixed",      function(d) stress_mgwrsar(d, H = stress_hdist(d))),
  stress_config("gwr_gauss_fixed_tiny", function(d) stress_mgwrsar(d, H = stress_hdist(d, 0.02))),
  stress_config("gwr_bisq_adapt30",     function(d) stress_mgwrsar(d, kernels = "bisq", H = min(30, nrow(d) - 1), adaptive = TRUE)),
  stress_config("gwr_bisq_adapt6",      function(d) stress_mgwrsar(d, kernels = "bisq", H = 6, adaptive = TRUE)),
  stress_config("gwr_gauss_fixed_SE",   function(d) stress_mgwrsar(d, H = stress_hdist(d), SE = TRUE)),
  stress_config("gwr_gauss_fixed_Shat", function(d) stress_mgwrsar(d, H = stress_hdist(d), get_s = TRUE)),
  stress_config("mgwr_gauss_fixed",     function(d) stress_mgwrsar(d, Model = "MGWR", H = stress_hdist(d), fixed_vars = "X3")),
  stress_config("mgwr_bisq_adapt30_SE", function(d) stress_mgwrsar(d, Model = "MGWR", kernels = "bisq", H = min(30, nrow(d) - 1), adaptive = TRUE, fixed_vars = "X3", SE = TRUE)),
  stress_config("gtwr_gauss_fixed",     function(d) stress_mgwrsar(d, kernels = c("gauss", "gauss"), H = c(stress_hdist(d), stress_htime(d)), Type = "GDT", adaptive = c(FALSE, FALSE))),
  stress_config("gtwr_gauss_fixed_tiny",function(d) stress_mgwrsar(d, kernels = c("gauss", "gauss"), H = c(stress_hdist(d, 0.02), stress_htime(d, 0.02)), Type = "GDT", adaptive = c(FALSE, FALSE))),
  stress_config("gtwr_bisq_adapt30",    function(d) stress_mgwrsar(d, kernels = c("bisq", "bisq"), H = c(min(30, nrow(d) - 1), min(30, nrow(d) - 1)), Type = "GDT", adaptive = c(TRUE, TRUE))),
  stress_config("gtwr_cyclic_gauss",    function(d) stress_mgwrsar(d, kernels = c("gauss", "gauss_SYM_365"), H = c(stress_hdist(d), 30), Type = "GDT", adaptive = c(FALSE, FALSE))),
  stress_config("twr_gauss_fixed",      function(d) stress_mgwrsar(d, kernels = "gauss", H = stress_htime(d), Type = "T")),
  stress_config("gtwr_hwb2010",         function(d) gtwr_HWB2010(formula = stress_formula, data = d, coords = stress_coords(d), time = d$time, h_st = stress_hdist(d) * sqrt(2), tau = (stress_hdist(d) / stress_htime(d))^2)),
  # --- bandwidth searches ---
  stress_config("search_gd_bisq_aicc",  function(d) search_bandwidths(formula = stress_formula, data = d, coords = stress_coords(d), kernels = "bisq", Model = "GWR",
                                            control = list(adaptive = TRUE, criterion = "AICc", Type = "GD", verbose = FALSE, ncore = 1), hs_range = c(6, nrow(d)), n_seq = 8, n_rounds = 2, ncore = 1, verbose = FALSE)$best_model, family = "search"),
  stress_config("search_gdt_gauss_aicc",function(d) search_bandwidths(formula = stress_formula, data = d, coords = stress_coords(d), kernels = c("gauss", "gauss"), Model = "GWR",
                                            control = list(Z = d$time, adaptive = c(FALSE, FALSE), criterion = "AICc", Type = "GDT", verbose = FALSE, ncore = 1),
                                            hs_range = c(stress_hdist(d, 0.01), stress_hdist(d, 0.99)), ht_range = c(stress_htime(d, 0.01), stress_htime(d, 0.99)), n_seq = 6, n_rounds = 2, ncore = 1, verbose = FALSE)$best_model, family = "search"),
  # --- top-down scale ---
  stress_config("tds_mgwr_gauss",       function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_mgwr", kernels = "gauss", fixed_vars = NULL,
                                            control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS"),
                                            control = list(NN = nrow(d), adaptive = FALSE, Type = "GD", ncore = 1)), family = "tds"),
  stress_config("tds_mgwr_bisq_adapt",  function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_mgwr", kernels = "bisq", fixed_vars = NULL,
                                            control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS"),
                                            control = list(NN = nrow(d), adaptive = TRUE, Type = "GD", ncore = 1)), family = "tds"),
  stress_config("tds_mgtwr_gauss",      function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_mgtwr", kernels = c("gauss", "gauss"), fixed_vars = NULL,
                                            control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS"),
                                            control = list(Z = d$time, NN = nrow(d), adaptive = c(FALSE, FALSE), Type = "GDT", ncore = 1)), family = "tds"),
  # --- prediction ---
  stress_config("predict_gwr_model",    function(d) stress_predict(d, function(a) stress_mgwrsar(a, H = stress_hdist(a)), method_pred = "model", beta_proj = TRUE, exposant = 8), family = "predict"),
  stress_config("predict_gwr_shepard",  function(d) stress_predict(d, function(a) stress_mgwrsar(a, kernels = "bisq", H = min(30, nrow(a) - 1), adaptive = TRUE), method_pred = "shepard"), family = "predict"),
  stress_config("predict_gtwr_model",   function(d) stress_predict(d, function(a) stress_mgwrsar(a, kernels = c("gauss", "gauss"), H = c(stress_hdist(a), stress_htime(a)), Type = "GDT", adaptive = c(FALSE, FALSE)), method_pred = "model", beta_proj = TRUE, exposant = 8), family = "predict"),
  stress_config("predict_tds_mgwr",     function(d) stress_predict(d, function(a) TDS_MGWR(formula = stress_formula, data = a, coords = stress_coords(a), Model = "tds_mgwr", kernels = "bisq", fixed_vars = NULL,
                                            control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS"), control = list(NN = nrow(a), adaptive = TRUE, Type = "GD", ncore = 1)),
                                            method_pred = "model", beta_proj = TRUE, exposant = 8), family = "predict")
)

# ---- phase 2: other estimators, post-processing, deliberately bad inputs ----

stress_W <- function(d, k = 8) kernel_matW(H = k, kernels = "rectangle", coords = stress_coords(d), NN = min(k + 2, nrow(d)),
                                           TP = seq_len(nrow(d)), Type = "GD", adaptive = TRUE, diagnull = TRUE)

# a "clean error": an explicit message, not an index / NA / dimension error
stress_expect_error <- function(expr, pattern) {
  res <- tryCatch({ force(expr); "no error" }, error = function(e) conditionMessage(e))
  inv <- character()
  if (identical(res, "no error")) inv <- "no error raised"
  else if (!grepl(pattern, res, ignore.case = TRUE)) inv <- paste("unclear error:", substr(res, 1, 100))
  list(invariants = inv)
}

stress_configs2 <- function() list(
  # --- other estimators ---
  stress_config("sar_2sls",             function(d) { W <- stress_W(d); MGWRSAR(formula = stress_formula, data = d, coords = stress_coords(d), fixed_vars = NULL, kernels = "gauss", H = 1, Model = "SAR", control = list(W = W, Method = "2SLS")) }, family = "fit2"),
  stress_config("mgwrsar_1_0_kv",       function(d) { W <- stress_W(d); MGWRSAR(formula = stress_formula, data = d, coords = stress_coords(d), fixed_vars = NULL, kernels = "bisq", H = min(30, nrow(d) - 3), Model = "MGWRSAR_1_0_kv", control = list(W = W, adaptive = TRUE)) }, family = "fit2"),
  stress_config("mgwrsar_0_kc_kv",      function(d) { W <- stress_W(d); MGWRSAR(formula = stress_formula, data = d, coords = stress_coords(d), fixed_vars = "X3", kernels = "bisq", H = min(30, nrow(d) - 3), Model = "MGWRSAR_0_kc_kv", control = list(W = W, adaptive = TRUE)) }, family = "fit2"),
  stress_config("gwr_glm_binomial",     function(d) { d$Yb <- as.numeric(d$Y > median(d$Y)); MGWRSAR(formula = Yb ~ X1 + X2 + X3, data = d, coords = stress_coords(d), fixed_vars = NULL, kernels = "bisq", H = min(40, nrow(d) - 3), Model = "GWR_glm", control = list(adaptive = TRUE, family = binomial())) }, family = "fit2"),
  stress_config("gtwr_hwb_tau0_tauInf", function(d) { m0 <- gtwr_HWB2010(formula = stress_formula, data = d, coords = stress_coords(d), time = d$time, h_st = stress_hdist(d) * sqrt(2), tau = 0); mi <- gtwr_HWB2010(formula = stress_formula, data = d, coords = stress_coords(d), time = d$time, h_st = stress_htime(d), tau = Inf); inv <- c(stress_invariants(m0, d), stress_invariants(mi, d)); list(invariants = inv) }, family = "fit2"),
  stress_config("multiscale_gwr_bisq",  function(d) multiscale_gwr(formula = stress_formula, data = d, coords = stress_coords(d), kernels = "bisq", control_mgwr = list(maxiter = 5, verbose = FALSE), control = list(adaptive = TRUE, NN = nrow(d), ncore = 1)), family = "fit2"),
  stress_config("atds_gwr_gauss",       function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "atds_gwr", kernels = "gauss", fixed_vars = NULL, control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS", nrounds = 2), control = list(NN = nrow(d), adaptive = FALSE, Type = "GD", ncore = 1)), family = "fit2"),
  stress_config("tds_mgwr_fixed_vars",  function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_mgwr", kernels = "bisq", fixed_vars = "X3", control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS"), control = list(NN = nrow(d), adaptive = TRUE, Type = "GD", ncore = 1)), family = "fit2"),
  stress_config("tds_mgwr_init_gwr",    function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_mgwr", kernels = "bisq", fixed_vars = NULL, control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "GWR"), control = list(NN = nrow(d), adaptive = TRUE, Type = "GD", ncore = 1)), family = "fit2"),
  stress_config("tds_mgtwr_cyclic",     function(d) TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_mgtwr", kernels = c("gauss", "gauss_SYM_365"), fixed_vars = NULL, control_tds = list(nns = 10, verbose = FALSE, tol = 1e-3, init_model = "OLS"), control = list(Z = d$time, NN = nrow(d), adaptive = c(FALSE, FALSE), Type = "GDT", ncore = 1)), family = "fit2"),
  # --- post-processing ---
  stress_config("summary_print_gwr",    function(d) { m <- stress_mgwrsar(d, kernels = "bisq", H = min(30, nrow(d) - 3), adaptive = TRUE, SE = TRUE); invisible(capture.output(print(summary(m)))); invisible(capture.output(print(m))); list(invariants = character()) }, family = "post"),
  stress_config("summary_mgwr_fdr",     function(d) { m <- stress_mgwrsar(d, Model = "MGWR", kernels = "bisq", H = min(30, nrow(d) - 3), adaptive = TRUE, fixed_vars = "X3", SE = TRUE); invisible(capture.output(print(summary(m, fdr_method = "spatial_BY")))); list(invariants = character()) }, family = "post"),
  stress_config("coef_residuals_fitted",function(d) { m <- stress_mgwrsar(d, H = stress_hdist(d)); inv <- character(); if (!isTRUE(all.equal(as.numeric(fitted(m)), as.numeric(m@fit)))) inv <- c(inv, "fitted() != @fit"); if (!isTRUE(all.equal(as.numeric(residuals(m)), as.numeric(m@residuals)))) inv <- c(inv, "residuals() != @residuals"); cf <- coef(m); if (is.null(cf)) inv <- c(inv, "coef() NULL"); list(invariants = inv) }, family = "post"),
  stress_config("predict_gwr_TP",       function(d) stress_predict(d, function(a) stress_mgwrsar(a, kernels = "bisq", H = min(30, nrow(a) - 3), adaptive = TRUE), method_pred = "TP"), family = "post"),
  stress_config("predict_gwr_tWtp",     function(d) stress_predict(d, function(a) stress_mgwrsar(a, H = stress_hdist(a)), method_pred = "tWtp_model"), family = "post"),
  stress_config("predict_mgwr_model",   function(d) stress_predict(d, function(a) stress_mgwrsar(a, Model = "MGWR", H = stress_hdist(a), fixed_vars = "X3"), method_pred = "model"), family = "post"),
  stress_config("predict_far_outside",  function(d) { n <- nrow(d); app <- 1:(n - 10); m <- stress_mgwrsar(d[app, ], H = stress_hdist(d[app, ])); nc <- stress_coords(d)[(n - 9):n, ] + 1e3 * diff(range(d$x)); colnames(nc) <- c("X", "Y"); p <- predict(m, newdata = d[(n - 9):n, ], newdata_coords = nc, method_pred = "model", beta_proj = TRUE, exposant = 8); y <- if (is.list(p)) p$Y_predicted else p; list(invariants = if (all(is.finite(y))) character() else "prediction far outside not finite") }, family = "post"),
  stress_config("bootstrap_test_B4",    function(d) { m0 <- stress_mgwrsar(d, kernels = "bisq", H = min(30, nrow(d) - 3), adaptive = TRUE); m1 <- stress_mgwrsar(d, Model = "MGWR", kernels = "bisq", H = min(30, nrow(d) - 3), adaptive = TRUE, fixed_vars = "X3"); r <- mgwrsar_bootstrap_test(m1, m0, B = 4, ncore = 1); list(invariants = if (is.list(r) || is.numeric(r)) character() else "unexpected bootstrap result") }, family = "post"),
  # --- deliberately bad inputs: the error must be explicit ---
  stress_config("bad_model_name",       function(d) stress_expect_error(stress_mgwrsar(d, Model = "GWRX", H = stress_hdist(d)), "Model"), family = "bad"),
  stress_config("bad_type",             function(d) stress_expect_error(stress_mgwrsar(d, H = stress_hdist(d), Type = "GDX"), "Type"), family = "bad"),
  stress_config("bad_kernel_name",      function(d) stress_expect_error(stress_mgwrsar(d, kernels = "gaus", H = stress_hdist(d)), "kernel"), family = "bad"),
  stress_config("bad_H_negative",       function(d) stress_expect_error(stress_mgwrsar(d, H = -1), "H|bandwidth"), family = "bad"),
  stress_config("bad_H_NA",             function(d) stress_expect_error(stress_mgwrsar(d, H = NA), "H|bandwidth"), family = "bad"),
  stress_config("bad_H_adapt_gt_n",     function(d) stress_expect_error(stress_mgwrsar(d, kernels = "bisq", H = nrow(d) + 5, adaptive = TRUE), "neighbours|sample"), family = "bad"),
  stress_config("bad_NN_lt_p",          function(d) { w <- character(); withCallingHandlers(stress_mgwrsar(d, kernels = "bisq", H = 2, adaptive = TRUE, extra = list(NN = 2)), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }); list(invariants = if (any(grepl("global OLS", w))) character() else "no warning about the global OLS fallback") }, family = "bad"),
  stress_config("bad_W_dims",           function(d) stress_expect_error(MGWRSAR(formula = stress_formula, data = d, coords = stress_coords(d), fixed_vars = NULL, kernels = "gauss", H = 1, Model = "SAR", control = list(W = stress_W(d[1:10, ]), Method = "2SLS")), "W|dimension"), family = "bad"),
  stress_config("bad_fixed_vars_name",  function(d) stress_expect_error(stress_mgwrsar(d, Model = "MGWR", H = stress_hdist(d), fixed_vars = "X9"), "fixed_vars|X9"), family = "bad"),
  stress_config("bad_gdt_without_Z",    function(d) stress_expect_error(MGWRSAR(formula = stress_formula, data = d, coords = stress_coords(d), fixed_vars = NULL, kernels = c("gauss", "gauss"), H = c(stress_hdist(d), 10), Model = "GWR", control = list(Type = "GDT", adaptive = c(FALSE, FALSE))), "Z|time"), family = "bad"),
  stress_config("bad_na_in_data",       function(d) { d$X1[3] <- NA; stress_expect_error(stress_mgwrsar(d, H = stress_hdist(d)), "NA|missing") }, family = "bad"),
  stress_config("bad_na_in_coords",     function(d) { d$x[3] <- NA; stress_expect_error(stress_mgwrsar(d, H = stress_hdist(d)), "NA|coords|missing") }, family = "bad"),
  stress_config("bad_tds_model",        function(d) stress_expect_error(TDS_MGWR(formula = stress_formula, data = d, coords = stress_coords(d), Model = "tds_xxx", kernels = "gauss", control_tds = list(nns = 5), control = list(adaptive = TRUE)), "Model"), family = "bad"),
  stress_config("bad_search_criterion", function(d) stress_expect_error(search_bandwidths(formula = stress_formula, data = d, coords = stress_coords(d), kernels = "bisq", Model = "GWR", control = list(adaptive = TRUE, criterion = "AICcc", Type = "GD", verbose = FALSE, ncore = 1), hs_range = c(6, nrow(d)), n_seq = 4, n_rounds = 1, ncore = 1, verbose = FALSE), "criterion"), family = "bad")
)
