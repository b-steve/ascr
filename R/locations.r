locations <- function(fit, id){
    traps <- fit$traps
    mask <- fit$mask
    n.mask <- nrow(mask)
    a <- attr(mask, "area")
    detfn <- fit$detfn
    det.pars <- getpar(fit, fit$detpars, as.list = TRUE)
    dists <- distances(traps, mask)
    ## Calculating density due to animal locations.
    p.det <- p.dot(fit)
    esa <- a*sum(p.dot(points = mask, traps = traps, detfn = detfn, pars = det.pars))
    f.x <- p.det/esa
    ## Calculating conditional density of capture history, given location.
    capt <- fit$capt$bincapt[id, ]
    det.probs <- calc.detfn(dists, detfn, det.pars)
    f.capt <- aaply(det.probs*capt + (1 - det.probs)*(1 - capt), 2, prod)
    ## Putting into a matrix for contours().
    unique.x <- sort(unique(mask[, 1]))
    unique.y <- sort(unique(mask[, 2]))
    z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
    for (i in 1:n.mask){
        x <- mask[, 1][i]
        y <- mask[, 2][i]
        index.x <- which(x == unique.x)
        index.y <- which(y == unique.y)
        z[index.x, index.y] <- f.x[i]*f.capt[i]
    }
    contour(unique.x, unique.y, z)
    text(traps, col = "red", labels = 1:nrow(traps))
}
