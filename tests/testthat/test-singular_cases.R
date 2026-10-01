# Nearly singular local fits (tools/singular_cases.R): leverages must stay in
# [0, 1] whatever the conditioning of the local systems, and leverages, fitted
# values and (where they are reproducible) coefficients must match the
# reference registry tests/testthat/_singular_hashes.csv, computed with
# mgwrsar 1.3.2 (tools/write_singular_hashes.R). The Q-free leverage path of
# 1.4 failed the first two cases (leverages up to 22).

library(testthat)
library(mgwrsar)
source('../../tools/singular_cases.R')
source('../../tools/check_hash_against_registry.R')

sing_registry <- read.csv('_singular_hashes.csv', stringsAsFactors = FALSE)
sing_cases <- singular_cases()

for (nm in names(sing_cases)) {
  test_that(paste("nearly singular fit:", nm), {
    stopifnot(requireNamespace("digest", quietly = TRUE))
    ref <- sing_registry[sing_registry$case == nm, ]
    expect_equal(nrow(ref), 1L)
    m <- sing_cases[[nm]]()
    s <- singular_summary(m)
    n <- s$n

    # invariants of a hat matrix diagonal
    expect_true(all(is.finite(s$TS)))
    expect_true(all(s$TS >= -1e-8 & s$TS <= 1 + 1e-8))
    expect_lte(s$tS, n)
    expect_equal(s$tS, sum(s$TS), tolerance = 1e-8)
    # a saturated fit gets an infinite AICc, the others a finite one
    if (n - 1 - s$tS <= 0) expect_identical(s$AICc, Inf) else expect_true(is.finite(s$AICc))
    # the hat matrix, when formed, carries the same diagonal
    if (!is.null(s$Shat_diag)) expect_equal(unname(s$Shat_diag), s$TS, tolerance = 1e-10)

    # agreement with the reference version
    expect_identical(hash_coef_matrix(matrix(s$TS), digits = ref$digits_TS), ref$hash_TS)
    expect_identical(hash_coef_matrix(matrix(s$fit), digits = ref$digits_fit), ref$hash_fit)
    if (!is.na(ref$digits_Betav))
      expect_identical(hash_coef_matrix(unname(s$Betav), digits = ref$digits_Betav), ref$hash_Betav)
    expect_equal(signif(s$tS, 8), ref$tS, tolerance = 1e-7)
    if (is.finite(ref$AICc)) expect_equal(signif(s$AICc, 8), ref$AICc, tolerance = 1e-7)
  })
}
