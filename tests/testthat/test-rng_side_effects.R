# Fitting functions seed their own generator for internal draws and must leave
# the user's generator (kind and state) exactly as they found it.

source('../../tools/singular_cases.R')

rng_untouched <- function(expr) {
  old_kind <- RNGkind()
  set.seed(42); before <- runif(3)
  set.seed(42); force(expr); after <- runif(3)
  expect_identical(RNGkind(), old_kind)
  expect_identical(after, before)
}

test_that("MGWRSAR (GD and GDT) and simu_multiscale leave the user's RNG untouched", {
  data(mydata)
  coords <- as.matrix(mydata[, c("x", "y")])
  rng_untouched(MGWRSAR(formula = "Y_gwr ~ X1 + X2 + X3", data = mydata, coords = coords,
                        fixed_vars = NULL, kernels = "gauss", H = 20, Model = "GWR",
                        control = list(SE = FALSE, adaptive = TRUE)))
  # a space-time fit: prep_d() sub-samples the distances with its own seed
  rng_untouched(singular_cases()$gdt_gauss_tiny())
  rng_untouched(simu_multiscale(n = 100, myseed = 3))
  # and the seed is also restored when no .Random.seed existed yet
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  simu_multiscale(n = 50, myseed = 3)
  expect_false(exists(".Random.seed", envir = globalenv()))
})

test_that("search_bandwidths leaves the caller's connections open", {
  data(mydata)
  coords <- as.matrix(mydata[1:300, c("x", "y")])
  tf <- tempfile(); con <- file(tf, open = "w"); on.exit({ try(close(con), silent = TRUE); unlink(tf) }, add = TRUE)
  writeLines("before", con)
  out <- capture.output(s <- search_bandwidths(formula = Y_gwr ~ X1 + X2 + X3, data = mydata[1:300, ], coords = coords,
    kernels = "bisq", Model = "GWR", control = list(adaptive = TRUE, criterion = "AICc", Type = "GD", verbose = FALSE, ncore = 1),
    hs_range = c(10, 200), n_seq = 4, n_rounds = 1, ncore = 1, verbose = FALSE))
  expect_true(isOpen(con))
  writeLines("after", con); close(con)
  expect_identical(readLines(tf), c("before", "after"))
})
