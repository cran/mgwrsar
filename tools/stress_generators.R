# Adversarial data generators for the stress campaign (tools/stress_cases.R,
# tests/testthat/test-stress_*.R). Each returns a data.frame with coordinates
# x, y, a time column `time` (days), a response Y and predictors X1..X3 (plus
# extra columns named in the generator), built so that one empirical condition
# known or suspected to hurt local estimators is present. All draws use a
# pinned Mersenne-Twister generator (see with_mt_seed in tools/singular_cases.R).

if (!exists("with_mt_seed")) {
  f <- if (file.exists("tools/singular_cases.R")) "tools/singular_cases.R" else "../../tools/singular_cases.R"
  source(f)
}

# the report data use F1..F3: give them the campaign names
stress_sites <- function(seed = 123) {
  d <- singular_data_sites(seed = seed)
  names(d)[match(c("F1", "F2", "F3"), names(d))] <- c("X1", "X2", "X3")
  d
}

stress_base <- function(n = 300, seed = 1, spread = 1) with_mt_seed(seed, {
  d <- data.frame(x = runif(n) * spread, y = runif(n) * spread, time = runif(n, 1, 1000))
  d$X1 <- rnorm(n); d$X2 <- rnorm(n); d$X3 <- rnorm(n)
  b1 <- 2 + d$x / spread; b2 <- sin(2 * pi * d$y / spread); b3 <- 0.5
  d$Y <- 1 + b1 * d$X1 + b2 * d$X2 + b3 * d$X3 + rnorm(n, 0, 0.5)
  d
})

# Every generator runs under its own pinned seed (stress_generators()): the
# draws made after stress_base() (rpois, sample, rnorm, rt, runif) used the
# session generator and gave other data in every process (2026-10-01).
stress_generators <- function(seed_offset = 1000L) {
  g <- stress_generators_raw()
  out <- lapply(seq_along(g), function(i) { f <- g[[i]]; s <- seed_offset + i; function() with_mt_seed(s, f()) })
  names(out) <- names(g)
  out
}

stress_generators_raw <- function() list(
  # --- geometry of the design points ---
  repeated_sites     = function() stress_sites(123),                           # ~2.7 obs per site, 200 m jitter
  exact_duplicates   = function() { d <- stress_base(150, 2); d <- rbind(d, d); d[order(d$x), ] },
  collinear_points   = function() { d <- stress_base(300, 3); d$y <- 0.3 * d$x + 0.1; d },   # all points on a line
  two_clusters       = function() { d <- stress_base(300, 4); far <- d$x > 0.5; d$x[far] <- d$x[far] + 1000; d$y[far] <- d$y[far] + 1000; d },
  isolated_point     = function() { d <- stress_base(300, 5); d$x[1] <- 50; d$y[1] <- 50; d },
  anisotropic_scale  = function() { d <- stress_base(300, 6); d$x <- d$x * 1e6; d$y <- d$y * 1e-3; d },  # x in 1e6, y in 1e-3
  lonlat_degrees     = function() { d <- stress_base(300, 7); d$x <- -5 + 10 * d$x; d$y <- 42 + 8 * d$y; d },
  tiny_n             = function() stress_base(12, 8),                        # n close to p
  # --- design matrix ---
  constant_predictor = function() { d <- stress_base(300, 9); d$X3 <- 1; d },
  near_collinear     = function() { d <- stress_base(300, 10); d$X2 <- d$X1 + 1e-7 * rnorm(nrow(d)); d },
  exact_collinear    = function() { d <- stress_base(300, 11); d$X2 <- 2 * d$X1; d },
  local_constant     = function() { d <- stress_base(300, 12); d$X3 <- as.numeric(d$x > 0.5); d },  # dummy constant in each half
  binary_sparse      = function() { d <- stress_base(300, 13); d$X3 <- as.numeric(seq_len(nrow(d)) <= 5); d },  # 5 ones only
  huge_scale         = function() { d <- stress_base(300, 14); d$X1 <- d$X1 * 1e8; d$X2 <- d$X2 * 1e-8; d },
  integer_predictors = function() { d <- stress_base(300, 15); d$X1 <- rpois(nrow(d), 3); d$X2 <- sample(0:1, nrow(d), TRUE); d },
  # --- response and noise ---
  heavy_tails        = function() { d <- stress_base(300, 16); d$Y <- d$Y + rt(nrow(d), 1.5); d },
  outliers           = function() { d <- stress_base(300, 17); d$Y[1:3] <- d$Y[1:3] + 1e3; d },
  heteroskedastic    = function() { d <- stress_base(300, 18); d$Y <- d$Y + rnorm(nrow(d), 0, 5 * d$x); d },
  constant_response  = function() { d <- stress_base(300, 19); d$Y <- 3; d },
  perfect_fit        = function() { d <- stress_base(300, 20); d$Y <- 1 + 2 * d$X1 - d$X2; d },   # no noise, global linear
  binary_response    = function() { d <- stress_base(300, 21); d$Y <- as.numeric(d$Y > median(d$Y)); d },
  # --- time axis ---
  constant_time      = function() { d <- stress_base(300, 22); d$time <- 100; d },
  few_dates          = function() { d <- stress_base(300, 23); d$time <- sample(c(10, 500, 990), nrow(d), TRUE); d },
  time_ties_sorted   = function() { d <- stress_base(300, 24); d$time <- rep(1:30, each = 10); d },
  cyclic_boundary    = function() { d <- stress_base(300, 25); d$time <- ifelse(runif(nrow(d)) < 0.5, runif(nrow(d), 1, 10), runif(nrow(d), 355, 365)); d },
  # --- combinations ---
  sites_collinear_X  = function() { d <- stress_sites(5); d$X2 <- d$X1 + 1e-6 * with_mt_seed(99, rnorm(nrow(d))); d }
)
