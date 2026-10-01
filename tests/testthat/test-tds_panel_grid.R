# Panel data (T observations per site): the non-adaptive spatial grid of
# TDS_MGWR must be built on DISTINCT locations, so that its floor is the
# distance to the 2nd site and not a fraction of the observations.
test_that("non-adaptive spatial grid counts distinct sites on a balanced panel", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "set RUN_LONG_TESTS=1")
  set.seed(3)
  N <- 36L; T <- 5L
  u <- runif(N); v <- runif(N)      # random sites: continuous distance distribution
  U <- rep(u, T); V <- rep(v, T); tau <- rep(seq_len(T), each = N)
  b2 <- 84 * U * V * (1 - U) * (1 - V); b3 <- 3 * (U - V); b4 <- 4 * sin(6 * (U + V))
  X1 <- rnorm(N * T); X2 <- rnorm(N * T); X3 <- rnorm(N * T)
  eta <- 1 + tau / T + b2 * X1 + b3 * X2 + b4 * X3
  d <- data.frame(Y = eta + rnorm(N * T, sd = sd(eta) / 3), X1, X2, X3)
  # panel_floor = 1: this test checks the count-to-distance conversion, not the
  # panel extension of the grid below the first site (test-tds_levers.R)
  m <- TDS_MGWR(Y ~ X1 + X2 + X3, Model = "tds_mgtwr", data = d, coords = cbind(U, V),
                kernels = c("gauss", "gauss"),
                control_tds = list(nns = 8, get_AIC = FALSE, verbose = FALSE, ncore = 1,
                                   init_model = "OLS", panel_floor = 1),
                control = list(Z = tau, adaptive = c(FALSE, FALSE), NN = N * T, Type = "GDT"))
  # reference: site-level pairwise distances, quantile 2/(N-1) = the default floor minv = 2
  ds <- as.numeric(dist(cbind(u, v)))
  floor_sites <- quantile(ds, probs = 2 / (N - 1), names = FALSE)
  floor_obs <- quantile(ds, probs = 2 / (N * T - 1), names = FALSE)   # what an observation count would give
  Vd <- as.numeric(m@V)
  expect_true(all(Vd > 0))
  expect_equal(min(Vd), floor_sites, tolerance = 0.25)
  expect_gt(min(Vd), 1.5 * floor_obs)      # an observation count would put the floor much lower
  # the top of the grid is the maximum site distance, the bandwidths are on the grid
  expect_true(all(as.numeric(m@H) >= min(Vd) - 1e-8))
})

# control_tds$V_dist: a spatial grid given directly in distances (fixed kernel)
# is used as is and may go below the distance to the nearest site. Given alone
# (no V), it must not be picked up by `control_tds$V` through partial matching.
test_that("V_dist sets a distance grid below the nearest-site distance", {
  skip_if_not(nzchar(Sys.getenv("RUN_LONG_TESTS")), "set RUN_LONG_TESTS=1")
  set.seed(3)
  N <- 36L; T <- 5L
  u <- runif(N); v <- runif(N)
  U <- rep(u, T); V <- rep(v, T); tau <- rep(seq_len(T), each = N)
  b2 <- 84 * U * V * (1 - U) * (1 - V); b3 <- 3 * (U - V); b4 <- 4 * sin(6 * (U + V))
  X1 <- rnorm(N * T); X2 <- rnorm(N * T); X3 <- rnorm(N * T)
  eta <- 1 + tau / T + b2 * X1 + b3 * X2 + b4 * X3
  d <- data.frame(Y = eta + rnorm(N * T, sd = sd(eta) / 3), X1, X2, X3)
  ds <- as.numeric(dist(cbind(u, v)))
  grid <- exp(seq(log(max(ds)), log(min(ds) / 2), length.out = 9))
  m <- TDS_MGWR(Y ~ X1 + X2 + X3, Model = "tds_mgtwr", data = d, coords = cbind(U, V),
                kernels = c("gauss", "gauss"),
                control_tds = list(nns = 8, V_dist = grid, get_AIC = FALSE, verbose = FALSE,
                                   ncore = 1, init_model = "OLS"),
                control = list(Z = tau, adaptive = c(FALSE, FALSE), NN = N * T, Type = "GDT"))
  Vd <- sort(as.numeric(m@V), decreasing = TRUE)
  expect_equal(min(Vd), min(ds) / 2)                    # floor below the nearest-site distance
  expect_true(all(grid[grid < max(ds)] %in% Vd))        # grid used as given
  expect_true(all(as.numeric(m@H) >= min(Vd) - 1e-8))
})
