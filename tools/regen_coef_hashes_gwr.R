# tools/regen_coef_hashes_gwr.R
# ---------------------------------------------------------------------------
# Regenerate the `configs_gwr` rows of tests/testthat/_coef_hashes.csv so they
# match THIS machine's floating-point output, replicating EXACTLY the fit path
# and coefficient matrix used by tests/testthat/test-Shat_fitted_GWR_model.R.
#
# WHY: the coefficient hashes are SHA-256 of Betav rounded to 10 digits, which
# differs by ~1e-10 across BLAS/LAPACK builds and machines. The reference must
# therefore be (re)generated on the machine/CI where the tests run.
#
# USAGE (from the package root, on the canonical machine):
#   pkgbuild::clean_dll(".")            # ensure a clean, current C++ build
#   devtools::load_all(".", recompile = TRUE)
#   source("tools/regen_coef_hashes_gwr.R")
#   devtools::test(filter = "Shat_fitted")   # should be green
#
# It preserves the configs_mgwr rows untouched (their test is independent).
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("digest", quietly = TRUE))
source("tools/test_init.R")
source("tools/configs_list_estimation.R")
source("tools/check_hash_against_registry.R")   # provides hash_coef_matrix()

csv <- file.path("tests", "testthat", "_coef_hashes.csv")
reg <- read.csv(csv, stringsAsFactors = FALSE)

# ---- replicate test-Shat_fitted_GWR_model.R setup exactly --------------------
configs <- configs_gwr
cfg0    <- configs[[1]]
data    <- setup_test_data_full(n = cfg0$n, lambda = cfg0$lambda,
                                config_beta = cfg0$config_beta,
                                config_snr  = cfg0$config_snr)$GD
df     <- data$mydata
coords <- as.matrix(data$coords)
fml    <- data$formula
W      <- kernel_matW(H = 4, kernels = "rectangle", coords = coords,
                      NN = 5, adaptive = TRUE, diagnull = TRUE)

fit_one <- function(cfg) {
  MGWRSAR(formula = fml, data = df, fixed_vars = cfg$fixed_vars, coords = coords,
          Model = cfg$Model, kernels = cfg$kernels, H = cfg$H,
          control = list(adaptive = cfg$adaptive, Type = cfg$Type,
                         NN = cfg$NN, get_s = FALSE, W = W))
}

errs <- character(0); nupd <- 0L
for (i in seq_len(nrow(reg))) {
  id <- reg$config_id[i]
  if (!(id %in% names(configs))) next            # keep configs_mgwr rows untouched
  cfg <- configs[[id]]
  m <- tryCatch(fit_one(cfg), error = function(e) e)
  if (inherits(m, "error")) {
    errs <- c(errs, paste(id, "::", conditionMessage(m))); next
  }
  B <- m@Betav
  if (isTRUE(cfg$SE) && length(m@sev) > 0) B <- B + m@sev  # mirrors the test exactly
  reg$n[i]      <- nrow(B)
  reg$p[i]      <- ncol(B)
  reg$digits[i] <- 10
  reg$hash[i]   <- hash_coef_matrix(B, digits = 10)
  nupd <- nupd + 1L
}

if (length(errs)) {
  cat("ERRORS (rows left unchanged):\n"); cat(errs, sep = "\n"); cat("\n")
}
write.csv(reg, csv, row.names = FALSE)
cat(sprintf("Regenerated %d configs_gwr rows; kept %d other rows. Wrote %s\n",
            nupd, sum(!(reg$config_id %in% names(configs))), csv))
