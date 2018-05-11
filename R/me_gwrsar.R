#' me_gwrsar
#' to be documented
#' @usage me_gwrsar(model, data, W_hat = NULL, Wtrue = NULL, DGP, Model, verbose = TRUE)
#' @param model  to be documented
#' @param data  to be documented
#' @param W_hat  to be documented
#' @param Wtrue  to be documented
#' @param DGP  to be documented
#' @param Model  to be documented
#' @param verbose  to be documented
#' @keywords internal
#' @return to be documented
me_gwrsar <-
function(model, data, W_hat = NULL, Wtrue = NULL, DGP, Model, verbose = TRUE) {
    n <- nrow(data)
    names_X = attr(terms(model.frame(model$formula, data)), "term.labels")
    names_X = c("Intercept", names_X)
    data$Intercept = rep(1, n)
    k = length(names_X)
    BETA = matrix(0, ncol = k, nrow = n)
    for (i in 0:(k - 1)) {
        BETA[, i + 1] <- data[, paste("Beta", i, sep = "")]
    }
    BETAV = model$Betav
    BETAC = model$Betac
    kv = ncol(BETAV)
    kc = length(BETAC)
    if (!is.null(model$fixed_vars)) {
        index_XC = which(names_X %in% model$fixed_vars)
        index_XV = (1:k)[-index_XC]
    }
    else {
        index_XC = NULL
        index_XV = (1:k)
    }
    if (Model == "SAR")
        index_XC = 1:(kc - 1)
    if (Model == "OLS")
        index_XC = 1:kc
    if (Model == "MGWRSAR_1_kc_0") {
        index_XC = 1:kc
        index_XV = NULL
    }
    cat("\n######### DGP = ", DGP, " #########\n")
    cat("\n######### model_DGP = ", Model, " #########\n")
    if (!(DGP %in% c("OLS", "MWGR", "GWR"))) {
        iW <- solve(Diagonal(n, 1) - data$lambda * Wtrue)
    }
    else iW <- Diagonal(n, 1)
    DBETA = matrix(0, ncol = k, nrow = n)
    for (i in 1:k) {
        DBETA[, i] <- rowSums(as.matrix(iW * BETA[, i]))
    }
    ME = matrix(0, ncol = k, nrow = 1)
    if (Model %in% c("SAR", "MGWRSAR_0_0_kv", "MGWRSAR_0_kc_kv"))
        iW_hat <- solve(Diagonal(n, 1) - as.numeric(BETAC[length(BETAC)]) *
            W_hat)
    if (Model %in% c("MGWRSAR_1_0_kv", "MGWRSAR_1_kc_kv", "MGWRSAR_1_kc_0"))
        iW_hat <- solve(Diagonal(n, 1) - BETAV[, ncol(BETAV)] *
            W_hat)
    if (Model %in% c("OLS", "MGWR", "GWR"))
        iW_hat <- Diagonal(n, 1)
    if (Model == "MGWRSAR_1_0_kv")
        kv = kv - 1
    if (Model %in% c("MGWRSAR_0_0_kv", "GWR", "MGWRSAR_1_0_kv")) {
        for (i in 1:kv) ME[i] <- mean(DBETA[, i] - rowSums(as.matrix(iW_hat *
            BETAV[, i])))
    }
    if (Model %in% c("MGWR", "MGWRSAR_0_kc_kv", "MGWRSAR_1_kc_kv")) {
        for (i in index_XV) {
            j = which(i == index_XV)
            ME[i] <- mean(DBETA[, i] - rowSums(as.matrix(iW_hat *
                BETAV[, j])))
        }
        for (i in index_XC) {
            j = which(i == index_XC)
            ME[i] <- mean(DBETA[, i] - rowSums(as.matrix(iW_hat *
                BETAC[j])))
        }
    }
    if (Model %in% c("MGWRSAR_1_kc_0", "OLS", "SAR")) {
        for (i in index_XC) ME[i] <- mean(DBETA[, i] - rowSums(as.matrix(iW_hat *
            BETAC[i])))
    }
    ME
}
