# Nearly singular local fits: the cases behind tests/testthat/test-singular_cases.R
# and the reference registry tests/testthat/_singular_hashes.csv (written by
# tools/write_singular_hashes.R from mgwrsar 1.3.2).
#
# Each case returns a fitted mgwrsar object from MGWRSAR() with fixed
# bandwidths. The data are built so that some local weighted least-squares
# systems are nearly or exactly rank deficient: repeated locations, nearly
# collinear predictors, a predictor constant within a neighbourhood, exact
# duplicates, exactly determined neighbourhoods.

# The generators pin the random number generator kind: MGWRSAR() and the
# bandwidth searches set the session RNG to L'Ecuyer-CMRG (restored on exit
# since 1.4.1, not before), so a plain set.seed() would give different data
# after the first fit of a session with earlier versions.
with_mt_seed <- function(seed, expr) {
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}

# Repeated locations (about 2.7 observations per site, within 200 m), five
# years, three predictors that are noisy copies of one signal (as in a model
# aggregation). This is the data of the GTWR bug report of 2026-10-01.
singular_data_sites <- function(n_sites = 150, n = 400, seed = 123) with_mt_seed(seed, {
  sites <- data.frame(x = runif(n_sites, 5e5, 8e5), y = runif(n_sites, 68e5, 71e5))
  obs <- data.frame(site = sample(n_sites, n, TRUE), year = sample(1:5, n, TRUE),
                    day = sample(105:320, n, TRUE))
  d <- data.frame(x = sites$x[obs$site] + runif(n, -200, 200),
                  y = sites$y[obs$site] + runif(n, -200, 200),
                  time = (obs$year - 1) * 365 + obs$day - 104)
  signal <- sin(d$x / 6e4) + cos(d$y / 8e4) + d$time / 1600
  d$Y  <- signal + rnorm(n, 0, 0.5)
  d$F1 <- signal + rnorm(n, 0, 0.6)
  d$F2 <- signal + rnorm(n, 0, 0.8)
  d$F3 <- signal + rnorm(n, 0, 1.2)
  d
})

# Distinct locations on a square, independent predictors, plus a regional
# dummy (constant inside each of four quadrants, hence collinear with the
# intercept in any neighbourhood that stays inside a quadrant).
singular_data_regions <- function(n = 300, seed = 11) with_mt_seed(seed, {
  d <- data.frame(x = runif(n), y = runif(n))
  d$X1 <- rnorm(n); d$X2 <- rnorm(n)
  d$R  <- as.numeric(d$x > 0.5) + 2 * as.numeric(d$y > 0.5)   # 0..3, constant per quadrant
  d$Y  <- 1 + d$X1 * (1 + d$x) + d$X2 * d$y + 0.5 * d$R + rnorm(n, 0, 0.3)
  d
})

singular_fit <- function(d, formula, H, kernels, Type = "GD", adaptive = FALSE,
                         Model = "GWR", fixed_vars = NULL, get_s = FALSE) {
  coords <- as.matrix(d[, c("x", "y")])
  ctl <- list(Type = Type, adaptive = adaptive, criterion = "AICc", get_s = get_s)
  if (Type == "GDT") ctl$Z <- d$time
  MGWRSAR(formula = formula, data = d, coords = coords, fixed_vars = fixed_vars,
          kernels = kernels, H = H, Model = Model, control = ctl)
}

singular_cases <- function() {
  ds <- singular_data_sites()
  f3 <- Y ~ F1 + F2 + F3
  # exact duplicates of every observation (same coordinates, X and Y)
  dd <- rbind(ds, ds); dd <- dd[order(dd$x), ]; rownames(dd) <- NULL
  # two predictors that are collinear up to 1e-7
  dc <- ds; dc$F2 <- dc$F1 + 1e-7 * with_mt_seed(5, rnorm(nrow(dc)))
  dr <- singular_data_regions()
  list(
    # the report: space-time gaussian kernels, tiny bandwidths
    gdt_gauss_tiny      = function() singular_fit(ds, f3, c(3000, 95), c("gauss", "gauss"), Type = "GDT", adaptive = c(FALSE, FALSE)),
    # same data, spatial kernel only, even smaller bandwidth
    gd_gauss_tiny       = function() singular_fit(ds, f3, 1500, "gauss"),
    # exactly determined neighbourhoods: p + 2 neighbours, compact kernel
    gd_bisq_6nn         = function() singular_fit(ds, f3, 6, "bisq", adaptive = TRUE),
    # exact duplicates with an adaptive compact kernel (ties in the distances)
    gd_bisq_duplicates  = function() singular_fit(dd, f3, 10, "bisq", adaptive = TRUE),
    # nearly collinear predictors, moderate bandwidth
    gd_gauss_collinear  = function() singular_fit(dc, Y ~ F1 + F2 + F3, 40000, "gauss"),
    # locally constant predictor (regional dummy), adaptive gaussian kernel
    gd_gauss_local_const = function() singular_fit(dr, Y ~ X1 + X2 + R, 25, "gauss", adaptive = TRUE),
    # mixed model (one constant coefficient) on the report data: mixed engine
    mgwr_gdt_tiny       = function() singular_fit(ds, f3, c(3000, 95), c("gauss", "gauss"), Type = "GDT", adaptive = c(FALSE, FALSE), Model = "MGWR", fixed_vars = "F3"),
    # the report case with the hat matrix requested (Q-forming path)
    gdt_gauss_tiny_shat = function() singular_fit(ds, f3, c(3000, 95), c("gauss", "gauss"), Type = "GDT", adaptive = c(FALSE, FALSE), get_s = TRUE)
  )
}

# Quantities compared across versions
singular_summary <- function(m) {
  list(Betav = as.matrix(m@Betav), TS = as.numeric(m@TS), tS = m@tS, AICc = m@AICc,
       fit = as.numeric(m@fit), n = length(m@fit),
       Shat_diag = if (length(m@Shat)) diag(as.matrix(m@Shat)) else NULL)
}

# Rounding digits used for the registry hashes (see tools/write_singular_hashes.R).
singular_digits <- function() {
  list(
    gdt_gauss_tiny       = c(TS = 8, fit = 8, Betav = NA),
    gd_gauss_tiny        = c(TS = 8, fit = 8, Betav = NA),  # near-singular: coefficient attribution depends on the pivoting
    gd_bisq_6nn          = c(TS = 8, fit = 8, Betav = 7),
    gd_bisq_duplicates   = c(TS = 8, fit = 8, Betav = 7),
    gd_gauss_collinear   = c(TS = 5, fit = 4, Betav = NA),
    gd_gauss_local_const = c(TS = 8, fit = 8, Betav = 7),
    mgwr_gdt_tiny        = c(TS = 8, fit = 8, Betav = NA),   # 1.3.2 fails on this case: reference from 1.4.1
    gdt_gauss_tiny_shat  = c(TS = 8, fit = 8, Betav = NA)
  )
}
