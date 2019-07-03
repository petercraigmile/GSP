
GSP.sim <- function (n, mu, cov.fun, dists, theta)
  ## ======================================================================
  ## By Peter F. Craigmile, pfc@stat.osu.edu
  ##
  ## Simulate 'n' realizations of a Gaussian stochastic process of
  ## length 'nrow(dists)' with mean 'mu' and a covariance function
  ## 'cov.fun'.
  ## ======================================================================
{

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

        drop(Z + mu)
    }
}

