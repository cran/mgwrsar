#' SSR
#' to be documented
#' @usage SSR(model)
#' @param model  to be documented
#' @keywords internal
#' @return to be documented
SSR <-
function (model)
{
    term1 = 0
    term2 = 0
    if (!is.null(model$Betav))
        term1 <- rowSums(model$XV * model$Betav)
    if (!is.null(model$Betac))
        term2 <- model$XC %*% as.matrix(model$Betac)
    residuals <- model$Y - term1 - term2
    sum(residuals^2)
}
