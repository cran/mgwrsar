#' compute_ts
#' to be documented
#' @usage compute_ts(S,Snew,,EXACT=TRUE,foldsl=NULL,iter=NULL)
#' @param S to be documented
#' @param Snew to be documented
#' @param foldsl to be documented
#' @noRd
#' @return to be documented
compute_ts <- function(S, Snew, foldsl) {
  n <- nrow(S)
  if (length(foldsl) == 1) {
    # Cas exact : on met à jour toute la matrice
    St <- S + Snew - eigenMapMatMult(S, Snew)
    tS <- sum(diag(St))
  } else {
    # Approximation par blocs diagonaux
    St <- matrix(0, nrow = n, ncol = n)  # matrice vide à remplir bloc par bloc
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

