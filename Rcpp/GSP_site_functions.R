
Edist <- function (x, y=x) {

    if (!any(class(x)=="matrix")) {
        stop("'x' must be a matrix")
    }
    if (!any(class(y)=="matrix")) {
        stop("'y' must be a matrix")
    }
    if (ncol(x) != ncol(y)) {
        stop("'x' and 'y' must have the same number of columns")
    }
    
    cppEdist(x,y)
}



make.grid <- function (list.of.coords) {

    cbind(apply(expand.grid(list.of.coords), 2, function (x) x))
}
