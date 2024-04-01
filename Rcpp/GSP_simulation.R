
GSP.sim <- function (n, mu, cov.fun, dists, theta, X, beta)
  ## ======================================================================
  ## By Peter F. Craigmile, pfc@stat.osu.edu
  ##
  ## Simulate 'n' realizations of a Gaussian stochastic process of
  ## length 'nrow(dists)' with mean 'mu' and a covariance function
  ## 'cov.fun'.  If 'mu' is missing calculate the mean using
  ## mu = X %*% beta.
  ## ======================================================================
{
    ## If 'mu' is missing calculate the mean using the design matrix
    if (missing(mu)) {

        mu <- drop(X %*% beta)
    }

    ## calculate the covariance matrix
    Sigma <- cov.fun(dists, theta)
    
    ## 'p' is the dimension of the covariance matrix
    p <- ncol(Sigma)

    R <- chol(Sigma, pivot = TRUE)

    R <- R[, order(attr(R, "pivot"))]

    Z <- crossprod(matrix(rnorm(n * p), nrow = p, byrow = TRUE), R)    

    if (missing(mu)) {

        drop(Z)
    } else {

        drop(t(t(Z) + mu))
    }
}

