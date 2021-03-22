#' Loop for multivariate estimation under hazard ratio ordering (experimental,
#' will possibly me modified)
#' 
#' @return 
#' A list containing the estimated conditional CDFs, estimation diagnostics,
#' and the constraint matrices for estimation.
#' 
#' @keywords internal
multivHazardLoop <- function(X, thresholds, nThr, weights, cpY, pars) {
    constr <- compOrd(X)
    N <- nrow(X)
    constrMat <- matrix(nrow = N, ncol = N, FALSE)
    constrMat[constr$paths] <- TRUE
    cdf <- matrix(ncol = nThr - 1, nrow = N, 0)
    
    A <- trReduc(constr$paths, N)
    nConstr <- nrow(A)
    l <- rep(0, nConstr)
    A <- Matrix::sparseMatrix(i = rep(seq_len(nConstr), 2), j = as.vector(A), 
      x = rep(c(-1, 1), each = nConstr), dims = c(nConstr, N))
    P <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = weights)
    i <- 1
    I <- nThr - 1
    conv <- vector("logical", I)
    
    q <- - weights * sapply(cpY, FUN = function(x) mean(thresholds[i] < x))
    qp <- osqp::osqp(P = P, q = q, A = A, l = l, pars = pars)
    sol <- qp$Solve()
    cdf[, 1] <- pmin(1, pmax(0, sol$x))
    conv[1] <- identical(sol$info$status, "maximum iterations reached")   
    
     if (I > 1) {
        for (i in 2:I) {
          sel <- which(cdf[, i - 1] > 0)
          N <- length(sel)
          constrMatTmp <- which(constrMat[sel, sel], arr.ind = TRUE)
          A <- trReduc(constrMatTmp, N)
          nConstr <- nrow(A)
          A <- Matrix::sparseMatrix(i = rep(seq_len(nConstr), 2), j = as.vector(A), 
            x = rep(c(-1, 1), each = nConstr), dims = c(nConstr, N))
          l <- rep(0, nConstr)
          P <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = weights[sel] * cdf[sel, i - 1])
          q <-  -weights[sel] * sapply(cpY[sel], FUN = function(x) mean(thresholds[i] < x))
          qp <- osqp::osqp(P = P, q = q, A = A, l = l, pars = pars)
          qp$WarmStart(x = cdf[sel, i - 1])
          sol <- qp$Solve()
          cdf[sel, i] <- cdf[sel, i - 1] * pmin(1, pmax(0, sol$x))
          conv[i] <- identical(sol$info$status, "maximum iterations reached")
        }

     }
    cdf <- 1 - cdf
    diagnostic <- list(
       precision = ifelse(I > 1, abs(min(diff(t(cdf)))), 0),
       convergence = mean(conv)
    )
    
    return(list(cdf = cdf, diagnostic = diagnostic, constr = constr))
}