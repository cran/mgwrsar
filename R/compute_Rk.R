#' compute_Rk
#' to be documented
#' @usage compute_Rk(Rk,Sk, St, foldsl)
#' @param Rk to be documented
#' @param Sk to be documented
#' @param St to be documented
#' @noRd
#' @return to be documented
compute_Rk<-function(Rk,Sk, St, foldsl) {
  n <- nrow(Sk)
# Approximation par blocs diagonaux
    Rkk <- matrix(0, nrow = n, ncol = n)  # matrice vide à remplir bloc par bloc
    for (i in seq_along(foldsl)) {
      idx <- foldsl[[i]]
      Sk_i  <- Sk[idx, idx]
      St_i  <- St[idx, idx]
      Rk_i  <- Rk[idx, idx]
      Rk_i <- eigenMapMatMult(Sk_i,Rk_i)  + Sk_i- eigenMapMatMult(Sk_i,St_i)
      Rkk[idx, idx] <- Rk_i
    }
  return(Rkk)
}
