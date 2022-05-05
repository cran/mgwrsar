#' cell
#' to be documented
#' @usage cell(q, xylim, ...)
#' @param q to be documented
#' @param xylim to be documented
#' @param ... to be documented
#' @noRd
#' @return to be documented
cell<-function(q, xylim, ...) {
  if (is(q,"quadtree")) f <- cell.quadtree else f <- cell.quadtree.leaf
  f(q, xylim, ...)
}
