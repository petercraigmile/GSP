

Edist <- function (x, y=x) {

    if (class(x) != "matrix") {
        stop("'x' must be a matrix")
    }
    if (class(y) != "matrix") {
        stop("'y' must be a matrix")
    }
    if (ncol(x) != ncol(y)) {
        stop("'x' and 'y' must have the same number of columns")
    }
    
    m <- nrow(x)
    n <- nrow(y)

    D <- matrix(0, m, n)

    tx <- t(x)

    for (j in 1:n) {
        
        D[,j] <- sqrt(colSums((tx - y[j,])^2))
    }

    D
}


make.grid <- function (list.of.coords) {

    cbind(apply(expand.grid(list.of.coords), 2, function (x) x))
}
