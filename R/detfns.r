## Returns a detection function.
get.detfn <- function(detfn){
    switch(detfn, hn = calc.hn, hr = calc.hr, th = calc.th,
           lth = calc.lth, ss = calc.ss, log.ss = calc.log.ss)
}

calc.hn <- function(d, pars){
    if (!identical(sort(names(pars)), c("g0", "sigma"))){
        stop("Argument 'pars' must have named components 'g0' and 'sigma'.")
    }
    g0 <- pars$g0
    sigma <- pars$sigma
    g0*exp(-(d^2/(2*sigma^2)))
}

calc.hr <- function(d, pars){
    if (!identical(sort(names(pars)), c("g0", "sigma", "z"))){
        stop("Argument 'pars' must have named components 'g0', 'sigma', and 'z'.")
    }
    g0 <- pars$g0
    sigma <- pars$sigma
    z <- pars$z
    g0*(1 - exp(-((d/sigma)^-z)))
}

calc.th <- function(d, pars){
    if (!identical(sort(names(pars)), c("scale", "shape"))){
        stop("Argument 'pars' must have named components 'scale' and 'shape'.")
    }
    scale <- pars$scale
    shape <- pars$shape
    0.5 - 0.5*erf(d/scale - shape)
}

calc.lth <- function(d, pars){
    if (!identical(sort(names(pars)), c("scale", "shape.1", "shape.2"))){
        stop("Argument 'pars' must have named components 'scale', 'shape.1', and 'shape.2'.")
    }
    scale <- pars$scale
    shape.1 <- pars$shape.1
    shape.2 <- pars$shape.2
    0.5 - 0.5*erf(shape.1 - exp(shape.2 - scale*d))
}

calc.ss <- function(d, pars){
    if (!identical(sort(names(pars)), c("b0.ss", "b1.ss", "cutoff", "sigma.ss"))){
        stop("Argument 'pars' must have named components 'b0.ss', 'b1.ss', 'sigma.ss', and 'cutoff'.")
    }
    b0.ss <- pars$b0.ss
    b1.ss <- pars$b1.ss
    sigma.ss <- pars$sigma.ss
    cutoff <- pars$cutoff
    1 - pnorm(cutoff, mean = b0.ss - b1.ss*d, sd = sigma.ss)
}

calc.log.ss <- function(d, pars){
    if (!identical(sort(names(pars)), c("b0.ss", "b1.ss", "cutoff", "sigma.ss"))){
        stop("Argument 'pars' must have named components 'b0.ss', 'b1.ss', 'sigma.ss', and 'cutoff'.")
    }
    b0.ss <- pars$b0.ss
    b1.ss <- pars$b1.ss
    sigma.ss <- pars$sigma.ss
    cutoff <- pars$cutoff
    1 - pnorm(cutoff, mean = exp(b0.ss - b1.ss*d), sd = sigma.ss)
}
