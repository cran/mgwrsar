#' Sharpen and row-normalise a prediction weight matrix
#'
#' Raises the weights of each row to the power \code{exposant} and
#' row-normalises them, after dividing each row by its largest weight so
#' that the nearest training point weighs 1 before the power is taken.
#' The result equals \code{normW(W^exposant)} whenever that expression is
#' finite, but it never underflows: with \code{normW(W^exposant)}, a target
#' whose weights are all small (Gaussian kernel with a small bandwidth, far
#' target) gets a row of zeros after the power and a non-finite prediction
#' (a few MGWR predictions at n = 5,000 in the Monte Carlo experiments).
#' Rows with no positive weight (compact kernel, target outside every
#' support) are left at zero and counted in the attribute \code{empty_rows}.
#'
#' @param W A dense or sparse (Matrix) weight matrix, one row per target.
#' @param exposant Sharpening exponent (\eqn{\gamma \ge 1}).
#' @return The sharpened, row-normalised matrix (same storage as the input,
#'   sparse matrices as \code{dgCMatrix}), with the attribute
#'   \code{empty_rows} (number of rows without any positive weight).
#' @noRd
.sharpen_rows <- function(W, exposant) {
  if (inherits(W, "sparseMatrix")) {
    W <- methods::as(methods::as(W, "CsparseMatrix"), "generalMatrix")
    n_empty <- nrow(W)
    if (length(W@x)) {
      ri <- W@i + 1L
      mx <- tapply(W@x, ri, max)
      rmax <- numeric(nrow(W))
      rmax[as.integer(names(mx))] <- mx
      sc <- rmax[ri]
      ok <- sc > 0
      W@x[ok] <- (W@x[ok] / sc[ok])^exposant
      W@x[!ok] <- 0
      n_empty <- sum(rmax <= 0)
    }
  } else {
    W <- as.matrix(W)
    rmax <- apply(W, 1, max)
    ok <- is.finite(rmax) & rmax > 0
    W[ok, ] <- (W[ok, , drop = FALSE] / rmax[ok])^exposant
    W[!ok, ] <- 0
    n_empty <- sum(!ok)
  }
  W <- normW(W)
  attr(W, "empty_rows") <- n_empty
  W
}

#' Nearest-neighbour fallback for prediction rows without any positive weight
#'
#' A target so far from every training point that its kernel weights all
#' underflow (Gaussian kernel with a small bandwidth) or vanish (compact
#' kernel) gets a row of zeros, hence a zero coefficient. Such rows receive
#' weight 1 on the nearest training point in the metric of the component:
#' spatial distance in units of its spatial bandwidth (the distance to the
#' k-th nearest training point of the target for an adaptive bandwidth)
#' and, for a space-time component, temporal distance in units of its
#' temporal bandwidth; a sentinel axis is ignored.
#'
#' @param W Row-normalised weight matrix (targets x training points).
#' @param coords Coordinates, training rows \code{S} then targets \code{O};
#'   columns 1:2 spatial, column 3 temporal for a space-time model.
#' @param S,O Row indices of the training points and of the targets.
#' @param myH Bandwidth (pair) of the component.
#' @param adaptive Logical (vector): adaptive spatial (and temporal) bandwidth.
#' @param V1,Vt1 Largest grid values (the sentinels) on each axis.
#' @return \code{W} with the empty rows replaced, and the attribute
#'   \code{fallback_rows} (their number).
#' @noRd
.nearest_fallback <- function(W, coords, S, O, myH, adaptive, V1, Vt1 = NULL) {
  rs <- if (inherits(W, "sparseMatrix")) Matrix::rowSums(W) else rowSums(W)
  er <- which(!is.finite(rs) | rs <= 0)
  n_fb <- length(er)
  if (!n_fb) {
    attr(W, "fallback_rows") <- 0L
    return(W)
  }
  cs <- coords[S, , drop = FALSE]
  use_s <- myH[1] < V1
  use_t <- length(myH) > 1 && !is.null(Vt1) && myH[2] < Vt1 && ncol(coords) >= 3
  for (o in er) {
    ct <- coords[O[o], ]
    d2 <- 0
    if (use_s) {
      ds2 <- (cs[, 1] - ct[1])^2 + (cs[, 2] - ct[2])^2
      hs <- if (isTRUE(adaptive[1])) sqrt(sort(ds2)[min(max(1L, round(myH[1])), length(ds2))]) else myH[1]
      if (!is.finite(hs) || hs <= 0) hs <- 1
      d2 <- d2 + ds2 / hs^2
    }
    if (use_t) {
      ht <- myH[2]
      if (!is.finite(ht) || ht <= 0) ht <- 1
      d2 <- d2 + (cs[, 3] - ct[3])^2 / ht^2
    }
    if (!use_s && !use_t) d2 <- (cs[, 1] - ct[1])^2 + (cs[, 2] - ct[2])^2
    j <- which.min(d2)
    W[o, ] <- 0
    W[o, j] <- 1
  }
  # set after the sub-assignments, which drop the attributes of a Matrix
  attr(W, "fallback_rows") <- n_fb
  W
}
