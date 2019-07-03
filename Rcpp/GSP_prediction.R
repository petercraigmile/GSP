



GSP.pred <- function (theta, z, X=NULL, cov.fun, dists, beta=NULL,
                      dists.between, dists.pred, X.pred=NULL,
                      calc.var=TRUE, calc.cov=FALSE, 
                      without.nugget=FALSE) {
    
    Sigma <- cov.fun(dists, theta)
    
    omega <- theta    
    if (without.nugget) {
        
        omega[length(omega)] <- 0
    }

    BSigma <- cov.fun(dists.between, omega)
        
    if (is.null(X)) {

        inner <- solve(Sigma, z)
        pred  <- drop(crossprod(BSigma, inner))
        
    } else {

        inner <- solve(Sigma, z - X %*% beta)
        pred  <- drop(X.pred %*% beta + crossprod(BSigma, inner))
    }

    if (calc.cov) {
        
        PSigma   <- cov.fun(dists.pred, omega)
        pred.cov <- PSigma - crossprod(solve(Sigma, BSigma), BSigma)

        if (calc.var) {

            pred.var <- diag(pred.cov)
            
        } else {

            pred.var <- NA
        }
        
    } else {
        
        pred.cov <- NA
        
        if (calc.var) {

            pred.var <- cov.fun(0, omega) - colSums(solve(Sigma, BSigma) * BSigma)
            
        } else {

            pred.var <- NA
        }
    }
    
    list(pred     = pred,
         pred.var = pred.var,
         pred.cov = pred.cov)
}





