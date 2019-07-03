
## ======================================================================
## 'GSP' stands for Gaussian stationary process, defined by
## covariances through distances.
## ======================================================================



GSP.exp <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Exponential covariance at distance 'h' and parameter values 'theta'.
    ## (Matern with nu=1/2).
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## On the original scale (untransformed),
    ## theta = (partial sill, range parameter) or
    ## theta = (partial sill, range parameter, nugget).
    ## ======================================================================
    
    p <- length(theta)

    if (p<2 || p>3) {

        stop("The length of 'theta' must be 2 or 3")
    }

    if (missing(transform.fun)) {

        drop(cppGSPexp(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPexp(rbind(h), transform.fun(theta), eps))
    }
}



GSP.Diggle <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Diggle covariance at distance 'h' and parameter values 'theta'
    ## (Matern with nu=3/2).
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## theta = (partial sill, range parameter) or
    ## theta = (partial sill, range parameter, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<2 || p>3) {

        stop("The length of 'theta' must be 2 or 3")
    }

    if (!missing(transform.fun)) {
        
         theta <- transform.fun(theta)
    }
    
    if (missing(transform.fun)) {
        
        drop(cppGSPDiggle(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPDiggle(rbind(h), transform.fun(theta), eps))
    }
}



GSP.Gaussian <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Gaussian covariance at distance 'h' and parameter values 'theta'
    ## (Matern with nu=infinity).
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## theta = (partial sill, range parameter) or
    ## theta = (partial sill, range parameter, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<2 || p>3) {

        stop("The length of 'theta' must be 2 or 3")
    }

    if (!missing(transform.fun)) {
        
         theta <- transform.fun(theta)
    }

    if (missing(transform.fun)) {

        drop(cppGSPGaussian(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPGaussian(rbind(h), transform.fun(theta), eps))
    }
}



GSP.powerexp <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Power exponential covariance at distance 'h' and parameter values 'theta'.
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## On the original scale (untransformed),
    ## theta = (partial sill, range parameter, exponent) or
    ## theta = (partial sill, range parameter, exponent, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<3 || p>4) {

        stop("The length of 'theta' must be 3 or 4")
    }

    if (missing(transform.fun)) {
        
        drop(cppGSPpowerexp(rbind(h), theta, eps))
        
    } else {
        
        drop(cppGSPpowerexp(rbind(h), transform.fun(theta), eps))
    }
}





GSP.Matern <- function (h, theta, transform.fun, eps=1e-12) {
    ## ======================================================================
    ## Matern covariance at distance 'h' and parameter values 'theta'.
    ## If 'transform.fun' is given transform 'theta' using this function.
    ##
    ## On the original scale (untransformed),
    ## theta = (partial sill, range parameter, smoothness parameter) or
    ## theta = (partial sill, range parameter, smoothness parameter, nugget).
    ## ======================================================================

    p <- length(theta)

    if (p<3 || p>4) {

        stop("The length of 'theta' must be 3 or 4")
    }

    if (!missing(transform.fun)) {
        
         theta <- transform.fun(theta)
    }

    h.over.rp <- h / theta[2]
    nu        <- theta[3]   

    if (p==3) {
        
        Cov <- ifelse(h > eps,
                      theta[1] * 2.0^(1.0-nu) * h.over.rp^nu / gamma(nu) * besselK(h.over.rp, nu),
                      theta[1])
    } else {

        Cov <- ifelse(h > eps,
                      theta[1] * 2.0^(1.0-nu) * h.over.rp^nu / gamma(nu) * besselK(h.over.rp, nu),
                      theta[1] + theta[4])
    }

    return(Cov)
}


