#' Compute Updated Cumulative Smoothing Matrix and Trace
#'
#' Updates the cumulative smoothing matrix S and computes its trace (effective
#' number of parameters) using exact or block-diagonal approximation.
#'
#' @usage compute_ts(S, Snew, foldsl)
#' @param S Current cumulative smoothing matrix (n x n).
#' @param Snew New smoothing matrix for the current variable.
#' @param foldsl List of index vectors defining the blocks for approximation.
#'   If length 1, exact computation is used.
#' @noRd
#' @return A list with S (updated cumulative matrix) and tS (trace of S).
compute_ts <- function(S, Snew, foldsl) {
  n <- nrow(S)
  if (length(foldsl) == 1) {
    # Exact case: update the entire matrix
    St <- S + Snew - eigenMapMatMult(S, Snew)
    tS <- sum(diag(St))
  } else {
    # Block-diagonal approximation
    St <- matrix(0, nrow = n, ncol = n)  # empty matrix to fill block by block
    tS <- 0

    for (i in seq_along(foldsl)) {
      idx <- foldsl[[i]]
      S_i     <- S[idx, idx]
      Snew_i  <- Snew[idx, idx]
      St_i    <- S_i + Snew_i - eigenMapMatMult(S_i, Snew_i)

      St[idx, idx] <- St_i
      tS <- tS + sum(diag(St_i))
    }
  }

  return(list(S = St, tS = tS))
}

