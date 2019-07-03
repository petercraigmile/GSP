


GSP.m2l <- function (theta, z, cov.fun, dists, X=NULL, transform.fun) {
    ## if X is NULL then the mean is assumed to be zero.

    Sigma <- cov.fun(dists, theta, transform.fun)
    
    chol.Sigma <- tryCatch(chol(Sigma), error = function(e) e)

    if (inherits(chol.Sigma, "error")) {

        Inf

    } else {

        n   <- length(z)

        Lz <- backsolve(chol.Sigma, matrix(z, nrow=n), transpose = TRUE)
        
        if (!is.null(X)) {
            
            LX <- backsolve(chol.Sigma, X, transpose = TRUE)
            
            eps <- z - X %*% .lm.fit(LX, Lz)$coef
            
            Lz <- backsolve(chol.Sigma, matrix(eps, nrow=n), transpose = TRUE)
        }

        res <- 2*sum(log(diag(chol.Sigma))) + n * log(2 * pi) + colSums(Lz^2)

        return(res)
    }
}




    

GSP.beta.hat <- function (theta, z, cov.fun, dists, X=NULL) {
    ## if X is NULL then the mean is assumed to be zero.

    if (is.null(X)) {
        
        NA
        
    } else {

        Sigma <- cov.fun(dists, theta)

        chol.Sigma <- tryCatch(chol(Sigma), error = function(e) e)

        if (inherits(chol.Sigma, "error")) {
            
            Inf
            
        } else {

            LX <- backsolve(chol.Sigma, X, transpose = TRUE)
            Lz <- backsolve(chol.Sigma, matrix(z, nrow=length(z)), transpose = TRUE)

            .lm.fit(LX, Lz)$coef
        }
    }
}


GSP.beta.hat.cov <- function (theta, cov.fun, dists, X=NULL) {
    ## if X is NULL then the mean is assumed to be zero.

    if (is.null(X)) {
        
        NA
        
    } else {

        Sigma <- cov.fun(dists, theta)

        chol.Sigma <- tryCatch(chol(Sigma), error = function(e) e)

        if (inherits(chol.Sigma, "error")) {
            
            Inf
            
        } else {

            LX <- backsolve(chol.Sigma, X, transpose = TRUE)

            solve(crossprod(LX))
        }
    }
}



GSP.mle <- function (theta, z, X=NULL, cov.fun, dists, hessian=FALSE) {
    ## is cov.fun is missing using an independent normal model

    if (!missing(cov.fun)) {

        out <- optim(log(theta), GSP.m2l, method="BFGS",
                     z=z, cov.fun=cov.fun, dists=dists, X=X,
                     transform.fun=exp, hessian=hessian)

        theta.hat <- exp(out$par)

        if (!is.null(X)) {

            beta.hat <- GSP.beta.hat(theta.hat, z, cov.fun, dists, X)
        }
        
        m2l <- out$value

        neg.hessian <- 0.5 * out$hessian ## remember we calculation -2 * loglik!
        
    } else {

        ols <- lm.fit(X, z)

        theta.hat <- mean(ols$residuals^2)        

        beta.hat <- ols$coef

        m2l <- -2*sum(dnorm(z, drop(X %*% beta.hat), sqrt(theta.hat), log=TRUE))

        if (hessian) {
            
            neg.hessian <- length(z)/(2*theta.hat^2)
        } else {
            
            neg.hessian <- NULL
        }
    }

    if (is.null(X)) {

        list(theta.hat = theta.hat,
             AIC       = m2l+2*length(theta.hat),
             X         = X,
             Hessian   = neg.hessian)
        
    } else {

        list(theta.hat = theta.hat,
             beta.hat  = beta.hat,
             AIC       = m2l+2*(length(theta.hat)+ncol(X)),
             X         = X,
             Hessian   = neg.hessian)
    }
}



