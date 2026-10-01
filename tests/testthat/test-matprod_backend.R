test_that("matprod backends agree and options(mgwrsar.matprod) is honoured", {
  set.seed(1)
  A <- matrix(rnorm(60 * 40), 60, 40, dimnames = list(NULL, paste0("x", 1:40)))
  B <- matrix(rnorm(40 * 30), 40, 30)
  ref <- unname(A %*% B)

  old <- options(mgwrsar.matprod = "blas")
  on.exit(options(old), add = TRUE)
  r_blas <- mgwrsar:::.mgwrsar_matprod(A, B)

  options(mgwrsar.matprod = "eigen")
  r_eigen <- mgwrsar:::.mgwrsar_matprod(A, B)

  # same shape and no dimnames whatever the backend
  expect_identical(dim(r_blas), c(60L, 30L))
  expect_identical(dim(r_eigen), c(60L, 30L))
  expect_null(dimnames(r_blas))
  expect_null(dimnames(r_eigen))
  expect_equal(r_blas, ref, tolerance = 1e-12)
  expect_equal(r_eigen, ref, tolerance = 1e-12)

  # "blas" is exactly R's product, "auto" is one of the two backends
  expect_identical(r_blas, ref)
  options(mgwrsar.matprod = "auto")
  r_auto <- mgwrsar:::.mgwrsar_matprod(A, B)
  expect_type(mgwrsar:::.mgwrsar_blas_is_optimised(), "logical")
  expect_identical(r_auto,
                   if (mgwrsar:::.mgwrsar_blas_is_optimised()) r_blas else r_eigen)

  # the compiled product refuses non-conformable matrices and takes integers
  options(mgwrsar.matprod = "eigen")
  expect_error(mgwrsar:::.mgwrsar_matprod(A, A), "non-conformable")
  expect_equal(mgwrsar:::.mgwrsar_matprod(matrix(1:4, 2), matrix(1:4, 2)),
               matrix(c(7, 10, 15, 22), 2))
})
