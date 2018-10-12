## Functions that automatically generate starting values for parameters.

## Not lifted from the secr package.
## Grabs the average recapture distance, or something.
autosigma <- function(args){
    easy.out <- attr(args$mask, "buffer")/4
    if (!args$same.traplocs){
        capt <- args$capt
        traps <- args$traps
        bincapt <- capt$bincapt
        ave.rc.dist <- function(x){
            trap.ids <- which(x == 1)
            if (length(trap.ids) > 1){
                rc.locs <- traps[trap.ids, ]
                rc.dists <- distances(rc.locs, rc.locs)
                w <- length(rc.dists[rc.dists > 0])
                out <- mean(rc.dists[rc.dists > 0])
            } else {
                out <- NA
                w <- NA
            }
            c(out, w)
        }
        mean.dists <- apply(bincapt, 1, ave.rc.dist)
        mean.dists <- mean.dists[, !is.na(mean.dists[1, ]), drop = FALSE]
        out <- sum(mean.dists[1, ]*mean.dists[2, ])/sum(mean.dists[2, ])
        if (!is.finite(out)){
            out <- easy.out
            }
    } else {
        out <- easy.out
    }
    out
}

## Write own pdot function.
autoD <- function(args){
    detfn <- args$detfn
    mask <- args$mask
    A <- args$A
    traps <- args$traps
    sv <- args$sv
    detpar.names <- args$detpar.names
    ss.link <- args$ss.opts$ss.link
    survey.length <- args$survey.length
    pars <- sv[detpar.names]
    if (any(detpar.names == "sigma.b0.ss")){
        pars$sigma.b0.ss <- 0
    }
    cutoff <- args$ss.opts$cutoff
    if (!is.null(cutoff)){
        pars$cutoff <- cutoff
    }
    esa <- p.dot(points = mask, esa = TRUE, traps = traps, detfn = detfn,
                 ss.link = ss.link, pars = pars, n.quadpoints = 8)
    ## HT-like estimator for D is n/esa.
    log(nrow(args$capt$bincapt)/(esa*survey.length))
}

`autoD.(Intercept)` <- autoD

autog0 <- function(args){
    0.95
}

autoz <- function(args){
    1
}

autosigma.toa <- function(args){
    0.0025
}

autokappa <- function(args){
    10
}

autob0.ss <- function(args){
    spherical <- identity
    do.call(args$ss.opts$ss.link, list(max(args$capt$ss)))
}

autosigma.b0.ss <- function(args){
    spherical <- identity
    0.01*do.call(args$ss.opts$ss.link, list(max(args$capt$ss)))
}

autob1.ss <- function(args){
    buffer <- attr(args$mask, "buffer")
    cutoff <- args$ss.opts$cutoff
    max.ss <- max(args$capt$ss)
    out <- (max.ss - cutoff)/(buffer/2)
    if (args$ss.opts$ss.link == "log"){
        out <- out/max.ss
    }
    out
}

autob2.ss <- function(args){
    0.1
}

autosigma.ss <- function(args){
    ss <- c(args$capt$ss)
    sd(ss[ss >= args$ss.opts$cutoff])
}

autoalpha <- function(args){
    2
}

autoshape <- function(args){
    2
}

autoscale <- function(args){
    if (args$detfn == "th"){
        out <- autosigma(args)
    } else if (args$detfn == "lth") {
        out <- (log(autoshape.1(args)) - log(autoshape.1(args) - 0.5))/autosigma(args)
    } else {
        stop("Detection function not recognised.")
    }
    out
}

autoshape.1 <- function(args){
    ## Magic number based on clever maths.
    0.809017
}

autoshape.2 <- function(args){
    log(autoshape.1(args) + autoscale(args)*autosigma(args))
}
