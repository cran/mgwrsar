#' Ensure Uniqueness of Coordinates or Values in 1D, 2D, or 3D Structures
#'
#' This function enforces uniqueness in numeric data structures of dimension
#' one, two, or three by applying minimal perturbations to duplicated values.
#' It is primarily intended to avoid degeneracy issues (e.g., duplicated
#' coordinates or time stamps) in spatial and spatio-temporal modeling.
#'
#' Depending on the dimensionality of the input, uniqueness is enforced as follows:
#' \itemize{
#'   \item \strong{1D}: duplicated values are slightly perturbed to ensure uniqueness;
#'   \item \strong{2D}: duplicated \eqn{(x, y)} coordinate pairs are resolved by perturbing
#'         the second coordinate;
#'   \item \strong{3D}: duplicated \eqn{(x, y)} pairs and duplicated temporal coordinates
#'         are handled separately, ensuring uniqueness in both space and time.
#' }
#'
#' The perturbation magnitude is proportional to machine precision and the
#' scale of the input data, ensuring numerical stability while preserving the
#' original structure as closely as possible.
#'
#' @param M A numeric vector, matrix, or data frame with 1, 2, or 3 columns.
#'   For matrices or data frames, columns are interpreted as spatial
#'   (\code{x}, \code{y}) and optional temporal (\code{t}) coordinates.
#'
#' @return
#' If \code{M} is a vector, a numeric vector of the same length is returned.
#' If \code{M} is a matrix or data frame, a numeric matrix with the same
#' dimensions is returned, with standardized column names:
#' \code{"t"} (1D), \code{"x", "y"} (2D), or \code{"x", "y", "t"} (3D).
#'
#' @details
#' This function is designed as a low-level utility for preprocessing spatial
#' or spatio-temporal inputs prior to kernel-based local regression or distance
#' computations, where duplicated locations or timestamps may cause numerical
#' singularities.
#'
#' @export
make_unique_by_structure <- function(M) {

  eps <- max(abs(M), na.rm = TRUE) * .Machine$double.eps * 10
  ## ----------------------------------------------------
  ## 0. Accept vectors, matrices OR data.frames
  ## ----------------------------------------------------
  input_is_vector <- is.vector(M) && !is.list(M)

  if (input_is_vector) {
    Mmat <- matrix(M, ncol = 1)
  } else if (is.data.frame(M)) {
    Mmat <- as.matrix(M)
  } else if (is.matrix(M)) {
    Mmat <- M
  } else {
    stop("Input must be a vector, matrix, or data.frame with numeric columns.")
  }

  if (!is.numeric(Mmat))
    stop("All columns must be numeric.")

  ncolM <- ncol(Mmat)

  ## Helper: assign names only for matrices
  set_dimnames <- function(X) {
    if (ncol(X) == 1) {
      colnames(X) <- "t"
    } else if (ncol(X) == 2) {
      colnames(X) <- c("x", "y")
    } else if (ncol(X) == 3) {
      colnames(X) <- c("x", "y", "t")
    }
    X
  }

  ## ====================================================
  ## CASE 1 – ONE COLUMN (vector or matrix 1-col)
  ## ====================================================
  if (ncolM == 1) {

    dup <- duplicated(Mmat[,1]) | duplicated(Mmat[,1], fromLast = TRUE)

    if (any(dup)) {
      r <- ave(which(dup), Mmat[dup,1], FUN = seq_along)
      Mmat[dup, 1] <-
        Mmat[dup, 1] + (r - 1) * eps
    }

    ## ---- return type based on input ----
    if (input_is_vector) {
      return(as.numeric(Mmat[,1]))  ## return vector
    } else {
      return(set_dimnames(Mmat))    ## return matrix with colnames
    }
  }

  ## ====================================================
  ## CASE 2 – TWO COLUMNS: unique (x,y) pairs
  ## ====================================================
  if (ncolM == 2) {

    keys <- paste(Mmat[,1], Mmat[,2], sep = "_")
    dup <- duplicated(keys) | duplicated(keys, fromLast = TRUE)

    if (any(dup)) {
      r <- ave(which(dup), keys[dup], FUN = seq_along)
      Mmat[dup, 2] <- Mmat[dup, 2] + (r - 1) * eps
    }

    return(set_dimnames(Mmat))
  }

  ## ====================================================
  ## CASE 3 – THREE COLUMNS: unique (x,y) and unique t
  ## ====================================================
  if (ncolM == 3) {

    keys12 <- paste(Mmat[,1], Mmat[,2], sep = "_")
    dup12  <- duplicated(keys12) | duplicated(keys12, fromLast = TRUE)
    dup3   <- duplicated(Mmat[,3]) | duplicated(Mmat[,3], fromLast = TRUE)

    if (any(dup12)) {
      r12 <- ave(which(dup12), keys12[dup12], FUN = seq_along)
      Mmat[dup12, 2] <- Mmat[dup12, 2] + (r12 - 1) * eps
    }

    if (any(dup3)) {
      r3 <- ave(which(dup3), Mmat[dup3,3], FUN = seq_along)
      Mmat[dup3, 3] <- Mmat[dup3, 3] + (r3 - 1) * eps
    }

    return(set_dimnames(Mmat))
  }

  ## ====================================================
  ## OTHERWISE
  ## ====================================================
  stop("This function only supports 1D, 2D, or 3D structures.")
}
