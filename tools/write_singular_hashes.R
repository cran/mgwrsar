# tools/write_singular_hashes.R
# Reference registry for the nearly singular cases of tools/singular_cases.R:
# tests/testthat/_singular_hashes.csv. Run it with the reference version of
# the package first in .libPaths(), e.g.
#   MGWRSAR_LIB=/path/to/lib_1.3.2 Rscript tools/write_singular_hashes.R
#   MGWRSAR_LIB=/path/to/lib_1.4.1 Rscript tools/write_singular_hashes.R mgwr_gdt_tiny
# A case name as argument restricts the update to that case; the other rows
# are kept. The rounding digits of each hashed quantity come from
# singular_digits(): coefficients of nearly singular fits are not reproducible
# across implementations (relative differences up to 1e-4 between 1.3.2 and
# 1.4.1), whereas fitted values and leverages are; a quantity with NA digits
# is not hashed for that case.
if (nzchar(Sys.getenv("MGWRSAR_LIB"))) .libPaths(c(normalizePath(Sys.getenv("MGWRSAR_LIB")), .libPaths()))
suppressPackageStartupMessages(library(mgwrsar))
source("tools/singular_cases.R")
source("tools/check_hash_against_registry.R")
only <- commandArgs(TRUE)
registry_file <- "tests/testthat/_singular_hashes.csv"

rows <- list()
for (nm in names(singular_cases())) {
  if (length(only) && !(nm %in% only)) next
  m <- singular_cases()[[nm]]()
  s <- singular_summary(m)
  d <- singular_digits()[[nm]]
  rows[[nm]] <- data.frame(
    case = nm, source = as.character(packageVersion("mgwrsar")),
    n = s$n, p = ncol(s$Betav),
    digits_TS = d["TS"], digits_fit = d["fit"], digits_Betav = d["Betav"],
    hash_TS    = hash_coef_matrix(matrix(s$TS), digits = d["TS"]),
    hash_fit   = hash_coef_matrix(matrix(s$fit), digits = d["fit"]),
    hash_Betav = if (is.na(d["Betav"])) NA_character_ else hash_coef_matrix(unname(s$Betav), digits = d["Betav"]),
    tS = signif(s$tS, 8), AICc = signif(s$AICc, 8),
    stringsAsFactors = FALSE)
  cat(sprintf("%-22s %s  tS=%.3f  AICc=%.3f\n", nm, rows[[nm]]$source, s$tS, s$AICc))
}
new <- do.call(rbind, rows)
if (file.exists(registry_file)) {
  old <- read.csv(registry_file, stringsAsFactors = FALSE)
  old <- old[!(old$case %in% new$case), ]
  new <- rbind(old, new)
}
new <- new[order(match(new$case, names(singular_cases()))), ]
write.csv(new, registry_file, row.names = FALSE)
cat("written:", registry_file, "(", nrow(new), "rows )\n")
