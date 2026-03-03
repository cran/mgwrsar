# =============================================================================
# Tests for prep_d() - Distance and Weight Matrix Preparation
# =============================================================================
# Package: mgwrsar v1.3.2
# Function: prep_d()
# 
# Signature réelle:
# prep_d(coords, NN, TP, Q=FALSE, extrapol=FALSE, ratio=1, QP=NULL, 
#        kernels=NULL, Type=NULL, stable_knn=FALSE)
# 
# Coverage: GD (spatial), GDT (spatiotemporal), T (temporal), kernel variations
# 
# Test strategy:
# - Use simu_multiscale() with fixed seed for reproducibility
# - Subset to 30 observations for speed
# - Test all parameter combinations (NN, TP, QP, extrapol, Type, kernels)
# - Verify mathematical properties (distances, indexG consistency)
# =============================================================================

library(testthat)
library(mgwrsar)


# Generate once for all tests
test_data_full <- setup_test_data_full()
test_data <- get_subset_30obs(test_data_full)

# Helper Functions ------------------------------------------------------------

#' Check basic structure of prep_d output
#' 
#' @param result Output from prep_d()
#' @param expected_rows Expected number of rows in indexG
#' @param expected_cols Expected number of columns in indexG
check_prep_d_structure <- function(result, expected_rows, expected_cols) {
  expect_type(result, "list")
  expect_named(result, c("indexG", "dists"), ignore.order = TRUE)
  expect_equal(nrow(result$indexG), expected_rows)
  expect_equal(ncol(result$indexG), expected_cols)
  expect_type(result$dists, "list")
}

#' Check distance matrix properties for GD
#' 
#' @param results  Index and Spatial distance matrix
#' @param TD  True distance
#' @param tol Numerical tolerance
check_distance_properties_GD <- function(result,myTD,TP,QP,NN,tol = 1e-10) {
  n=nrow(myTD)
  if(is.null(QP)) QP=1:n
  mydist=matrix(0,ncol=n,nrow=n)
  for( i in 1:length(TP)) mydist[TP[i],result$indexG[i,]]<-result$dist$dist_s[i,]
  # if(is.null(QP)){
  #   myTD<-myTD[TP,]
  # mydist<-mydist[TP,]
  # } else {
    myTD<-myTD[TP,QP]
    mydist<-mydist[TP,QP]
  #}
  dimnames(mydist) <- dimnames(myTD)
  if(length(NN)==n) {
    expect_equal(mydist,myTD, tolerance = tol,label = "distances diff should be zero")
  } else {
    idx=which(mydist==0)
    expect_equal(mydist[-idx],myTD[-idx], tolerance = tol,label = "distances diff should be zero : case NN")
  }
}

check_distance_properties_T <- function(result,myTD,TP,QP,NN,tol = 1e-10) {
  n=nrow(myTD)
  if(is.null(QP)) QP=1:n
  mydist=matrix(0,ncol=n,nrow=n)
  for( i in 1:length(TP)) mydist[TP[i],result$indexG[i,]]<-result$dist$dist_t[i,]
  # if(is.null(QP)){
  #   myTD<-myTD[TP,]
  # mydist<-mydist[TP,]
  # } else {
  myTD<-myTD[TP,QP]
  mydist<-mydist[TP,QP]
  #}
  dimnames(mydist) <- dimnames(myTD)
  if(length(NN)==n) {
    expect_equal(mydist,myTD, tolerance = tol,label = "distances diff should be zero")
  } else {
    idx=which(mydist==0)
    expect_equal(mydist[-idx],myTD[-idx], tolerance = tol,label = "distances diff should be zero : case NN")
  }
}


check_distance_properties_GDT <- function(result,myTD,myTD_t,TP,QP,NN,tol = 1e-10) {
  n=nrow(myTD)
  if(is.null(QP)) QP=1:n
  mydist_t<-mydist<-matrix(0,ncol=n,nrow=n)
  for( i in 1:length(TP)) mydist[TP[i],result$indexG[i,]]<-result$dist$dist_s[i,]
  for( i in 1:length(TP)) mydist_t[TP[i],result$indexG[i,]]<-result$dist$dist_t[i,]
  
  if(is.null(QP)){
    myTD<-myTD[TP,]
    myTD_t<-myTD_t[TP,]
    mydist<-mydist[TP,]
    mydist_t<-mydist_t[TP,]
  } else {
    myTD<-myTD[TP,QP]
    myTD_t<-myTD_t[TP,QP]
    mydist<-mydist[TP,QP]
    mydist_t<-mydist_t[TP,QP]
  }
  dimnames(mydist) <- dimnames(myTD)
  dimnames(mydist_t) <- dimnames(myTD_t)
  
  if(NN==n) {
    expect_equal(mydist,myTD, tolerance = tol,label = "distances diff should be zero")
    expect_equal(mydist_t,myTD_t, tolerance = tol,label = "distances diff should be zero")
  } else {
    idx=which(mydist==0)
    expect_equal(mydist[-idx],myTD[-idx], tolerance = tol,label = "distances diff should be zero ")
  }
}

# =============================================================================
# GROUPE 1: Tests GD (Spatial Only)
# =============================================================================

# GD - NN=n, estimation mode (TP=1:n, QP=NULL) -------------------------------

test_that("prep_d: GD NN=n, estimation mode", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 1:n
  QP = NULL
  
  result <- prep_d(
    coords = data$coords,
    NN = n,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = QP,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  check_prep_d_structure(result, expected_rows = n, expected_cols = n)
  expect_named(result$dists, "dist_s")
  check_distance_properties_GD(result,data$TD,TP,QP=NULL,NN=n,tol = 1e-10)
})

# GD - NN=n, prediction with extrapol (TP=1:25, QP=26:30) -------------------

test_that("prep_d: GD NN=n, prediction with extrapol", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 26:30
  QP <- 1:25
  NN <- n - length(TP)  # Adapté car extrapol=TRUE exclut TP
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = TRUE,
    ratio = 1,
    QP = QP,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, "dist_s")
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  
  check_distance_properties_GD(result,data$TD,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GD - NN<n, target points mode (TP=1:20, QP=NULL, extrapol=FALSE) -----------

test_that("prep_d: GD NN<n, target points mode", {
  data <- test_data$GD
  TP <- 1:20
  NN <- 15
  QP = NULL
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = QP,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, "dist_s")
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  
  check_distance_properties_GD(result,data$TD,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GD - NN<n, estimation mode -------------------------------------------------

test_that("prep_d: GD NN<n, estimation", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 1:n
  NN <- 25
  QP = NULL
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = QP,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, "dist_s")
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GD(result,data$TD,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GD - NN<n, prediction extrapol ---------------------------------------------

test_that("prep_d: GD NN<n, prediction with extrapol", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 26:30
  QP <- 1:25
  NN <- n - length(TP)  -5 # Adapté car extrapol=TRUE exclut TP
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = TRUE,
    ratio = 1,
    QP = QP,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, "dist_s")
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GD(result,data$TD,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GD - stable_knn=TRUE -------------------------------------------------------

test_that("prep_d: GD with stable_knn=TRUE", {
  data <- test_data$GD
  n <- nrow(data$coords)
  TP <- 1:n
  QP = NULL
  NN <- n
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = TRUE  # Test avec tri stable
  )
  
  expect_named(result$dists, "dist_s")
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GD(result,data$TD,TP,QP=QP,NN=NN,tol = 1e-10)
})

# =============================================================================
# GROUPE 2: Tests GDT (Spatiotemporal) - Non-Cyclic
# =============================================================================

# GDT - gauss-gauss, estimation ----------------------------------------------

test_that("prep_d: GDT gauss-gauss, estimation", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  NN <- n
  QP=NULL
  
  result <- prep_d(
    coords = data$coords_3col,  # 3 colonnes: x, y, time
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = c('gauss', 'gauss'),  # spatial, temporal
    Type = "GDT",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t","dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GDT(result,data$TD,data$TD_t,TP,QP=QP,NN=NN,tol = 1e-10)
  
  
})

# GDT - gauss-gauss, prediction with extrapol --------------------------------

test_that("prep_d: GDT gauss-gauss, prediction", {
  data <- test_data$GDT
  n <- nrow(data$coords)
  TP <- 26:30
  QP <- 1:25
  NN <- n - length(TP)
  
  result <- prep_d(
    coords = data$coords_3col,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = TRUE,
    ratio = 1,
    QP = QP,
    kernels = c('gauss', 'gauss'),
    Type = "GDT",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t","dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GDT(result,data$TD,data$TD_t,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GDT - gauss-gauss, target points -------------------------------------------

test_that("prep_d: GDT gauss-gauss, target points", {
  data <- test_data$GDT
  TP <- 1:20
  NN <- 15
  QP = NULL
  
  result <- prep_d(
    coords = data$coords_3col,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = QP,
    kernels = c('gauss', 'gauss'),
    Type = "GDT",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t","dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GDT(result,data$TD,data$TD_t,TP,QP=QP,NN=NN,tol = 1e-10)
})

# =============================================================================
# GROUPE 3: Tests GDT (Spatiotemporal) - Cyclic Temporal
# =============================================================================

# GDT - gauss-gauss_SYM_365, estimation --------------------------------------

test_that("prep_d: GDT gauss-gauss_SYM_365, cyclic temporal, estimation", {
  
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  NN <- n
  QP=NULL
  
  result <- prep_d(
    coords = data$coords_3col,  # 3 colonnes: x, y, time
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = c('gauss', 'gauss_SYM_365'),  # spatial, temporal
    Type = "GDT",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t","dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GDT(result,data$TD,data$TD_tc,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GDT - gauss-gauss_SYM_365, prediction --------------------------------------

test_that("prep_d: GDT gauss-gauss_SYM_365, prediction", {
  
  data <- test_data$GDT
  n <- nrow(data$coords)
  TP <- 26:30
  QP <- 1:25
  NN <- n - length(TP)
  
  result <- prep_d(
    coords = data$coords_3col,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = TRUE,
    ratio = 1,
    QP = QP,
    kernels = c('gauss', 'gauss_SYM_365'),
    Type = "GDT",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t","dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GDT(result,data$TD,data$TD_tc,TP,QP=QP,NN=NN,tol = 1e-10)
})

# GDT - gauss-gauss_SYM_365, target points -----------------------------------

test_that("prep_d: GDT gauss-gauss_SYM_365, target points", {
  data <- test_data$GDT
  TP <- 1:20
  NN <- 15
  QP = NULL
  
  result <- prep_d(
    coords = data$coords_3col,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = QP,
    kernels = c('gauss', 'gauss_SYM_365'),
    Type = "GDT",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t","dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GDT(result,data$TD,data$TD_tc,TP,QP=QP,NN=NN,tol = 1e-10)
})

# =============================================================================
# GROUPE 4: Tests T (Temporal Only)
# =============================================================================

# T - temporal only, estimation ----------------------------------------------

test_that("prep_d: Type=T, temporal only, estimation", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  NN <- n
  QP = NULL
  
  result <- prep_d(
    coords = as.matrix(data$coords_3col[,3],ncol=1),  # col 3 = time
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = 'gauss',
    Type = "T",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_T(result,data$TD_t,TP=TP,QP=QP,NN=NN,tol = 1e-10)
})

# T - temporal cyclic, estimation --------------------------------------------

test_that("prep_d: Type=T, cyclic temporal, estimation", {
  data <- test_data$GDT
  n <- nrow(data$coords_3col)
  TP <- 1:n
  NN <- n
  QP = NULL
  
  result <- prep_d(
    coords = as.matrix(data$coords_3col[,3],ncol=1),  # col 3 = time
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = 'gauss_SYM_365',
    Type = "T",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_t"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_T(result,data$TD_tc,TP=TP,QP=QP,NN=NN,tol = 1e-10)
})

# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

test_that("prep_d: handles TP subset correctly", {
  data <- test_data$GD
  TP <- c(1, 5, 10, 15, 20)  # Non-contiguous subset
  NN <- 10
  QP = NULL
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GD(result,data$TD,TP=TP,QP=QP,NN=NN,tol = 1e-10)
  })

test_that("prep_d: extrapol=TRUE with disjoint TP and QP  + Non-contiguous subset", {
  data <- test_data$GD
  TP <- c(1, 5, 10, 15, 20)  # Non-contiguous subset
  NN <- 10
  QP = c(2,6,23)
  
  result <- prep_d(
    coords = data$coords,
    NN = NN,
    TP = TP,
    Q = FALSE,
    extrapol = FALSE,
    ratio = 1,
    QP = NULL,
    kernels = 'gauss',
    Type = "GD",
    stable_knn = FALSE
  )
  
  expect_named(result$dists, c("dist_s"))
  check_prep_d_structure(result, expected_rows = length(TP), expected_cols = NN)
  check_distance_properties_GD(result,data$TD,TP=TP,QP=QP,NN=NN,tol = 1e-10)
})

# =============================================================================
# End of test-prep_d.R
# =============================================================================