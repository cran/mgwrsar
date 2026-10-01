#' Compute Updated Rk Matrix for Backfitting Diagnostics
#'
#' Updates the Rk matrix used in the backfitting procedure for computing
#' effective degrees of freedom, using block-diagonal approximation.
#'
#' @usage compute_Rk(Rk, Sk, St, foldsl)
#' @param Rk Current Rk matrix (n x n) for a given variable.
#' @param Sk New smoothing matrix for the current variable.
#' @param St Cumulative smoothing matrix.
#' @param foldsl List of index vectors defining the blocks for approximation.
#' @noRd
#' @return Updated Rk matrix (or list with Rkk and vk for the variance version).
compute_Rk_with_Var <- function(Rk, Sk, St, foldsl) {
n <- nrow(Sk)
Rkk <- matrix(0, nrow = n, ncol = n)
vk <- numeric(n) # Vector to store rowSums(Rk^2)

for (i in seq_along(foldsl)) {
  idx <- foldsl[[i]]
  Sk_i <- Sk[idx, idx]
  St_i <- St[idx, idx]
  Rk_i <- Rk[idx, idx]

  # Update Rk (current formula)
  Rk_new_i <- .mgwrsar_matprod(Sk_i, Rk_i) + Sk_i - .mgwrsar_matprod(Sk_i, St_i)

  Rkk[idx, idx] <- Rk_new_i

  # CRUCIAL COMPUTATION: Sum of squared row values of the block
  # rowSums(Rk_new_i^2)
  vk[idx] <- rowSums(Rk_new_i^2)
}
return(list(Rkk = Rkk, vk = vk))
}
compute_Rk<-function(Rk,Sk, St, foldsl) {
  n <- nrow(Sk)
# Block-diagonal approximation
    Rkk <- matrix(0, nrow = n, ncol = n)  # empty matrix to fill block by block
    for (i in seq_along(foldsl)) {
      idx <- foldsl[[i]]
      Sk_i  <- Sk[idx, idx]
      St_i  <- St[idx, idx]
      Rk_i  <- Rk[idx, idx]
      Rk_i <- .mgwrsar_matprod(Sk_i,Rk_i)  + Sk_i- .mgwrsar_matprod(Sk_i,St_i)
      Rkk[idx, idx] <- Rk_i
    }
  return(Rkk)
}
