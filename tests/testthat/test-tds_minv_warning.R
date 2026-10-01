# control_tds$minv above the identifiability floor of a Gaussian kernel:
# TDS_MGWR must warn at setup, and again when a bandwidth stops on that floor.
test_that("minv above the Gaussian floor warns at setup and on a censored bandwidth", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "set RUN_LONG_TESTS=1")
  simu <- simu_multiscale(n = 300, myseed = 1, config_snr = 0.9,
                          config_beta = "spatiotemp", config_eps = "normal")
  d <- simu$mydata; co <- as.matrix(simu$coords)
  fit <- function(ct) TDS_MGWR(Y ~ X1 + X2 + X3, Model = "tds_mgwr", data = d, coords = co,
                               kernels = "gauss", control_tds = ct,
                               control = list(adaptive = FALSE, NN = 300, Type = "GD"))
  base <- list(nns = 10, get_AIC = FALSE, verbose = FALSE, init_model = "OLS")
  # default floor (2 neighbours): no warning
  expect_no_warning(m0 <- fit(base))
  # raised floor: setup warning, and the purely spatial coefficients (X2, X3)
  # end on the floor -> second warning naming them
  w <- character()
  m6 <- withCallingHandlers(fit(c(base, list(minv = 6))),
                            warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 2)
  expect_match(w[1], "above the identifiability floor of a Gaussian kernel")
  expect_match(w[2], "stopped on the grid floor")
  expect_true(all(c("X2", "X3") %in% names(m6@H)[m6@H <= min(m6@V)]))
  # the raised floor is indeed higher than the default one
  expect_gt(min(m6@V), min(m0@V))
})
