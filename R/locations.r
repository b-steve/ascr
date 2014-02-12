#' Plotting estimated locations
#'
#' Plots estimated densities of animal locations, which are latent
#' variables in SECR models.
#'
#' @param fit A fitted model from \link[admbsecr]{admbsecr}.
#' @param id A numeric vector with row numbers from \code{fit$capt},
#' indicating which individuals' locations are to be plotted.
#' @param infotypes A character vector indicating the type(s) of
#' information to be used when plotting the estimated density of
#' location.  Elements can be a subset of \code{"capt"}, \code{"ang"},
#' \code{"dist"}, \code{"toa"}, and \code{"combined"}, where
#' \code{"combined"} combines all information types together. When
#' signal strength information is used in the model fit, then
#' selecting \code{"capt"} here will use this information.
#' @param xlim A numeric vector of length 2, giving the x coordinate range.
#' @param ylim A numeric vector of length 2, giving the y coordinate range.
#' @param
locations <- function(fit, id, infotypes = "combined",
                      xlim = range(mask[, 1]),
                      ylim = range(mask[, 2]),
                      cols = list(combined = "black", capt = "purple",
                          ang = "green", dist = "red", toa = "blue"),
                      mask = fit$mask, add = FALSE){
    ## Setting up plotting area.
    if (!add){
        plot.new()
        plot.window(xlim = xlim, ylim = ylim, asp = 1)
        box()
        axis(1)
        axis(2)
    }
    ## Working out which contours to plot.
    plot.types <- c("combined", "capt", "ang", "dist", "toa") %in% infotypes
    names(plot.types) <- c("combined", "capt", "ang", "dist", "toa")
    ## Setting all to TRUE if combined is used.
    ## Some error catching.
    for (i in c("ang", "dist", "toa")){
        if (plot.types[i] & !fit$fit.types[i]){
            msg <- paste("Contours for information type '", i, "' cannot be plotted as this information was not used in the model 'fit'", sep = "")
            warning(msg)
            plot.types[i] <- FALSE
        }
    }
    n.mask <- nrow(mask)
    detfn <- fit$detfn
    det.pars <- getpar(fit, fit$detpars, as.list = TRUE)
    dists <- distances(fit$traps, mask)
    ## Calculating density due to animal locations.
    p.det <- p.dot(fit = fit, points = mask)
    ## Divide by normalising constant; not conversion to square metres.
    a <- attr(mask, "area")
    f.x <- p.det/(a*sum(p.det))
    ## Calculating conditional density of capture history, given location.
    for (i in id){
        if (plot.types["combined"]){
            f.combined <- f.x
        }
        ## Contour due to capture history.
        if (plot.types["capt"] | plot.types["combined"]){
            capt <- fit$capt$bincapt[i, ]
            det.probs <- calc.detfn(dists, detfn, det.pars)
            f.capt <- colProds(det.probs*capt + (1 - det.probs)*(1 - capt))
            if (plot.types["capt"]){
                show.contour(mask, f.x*f.capt, cols$capt)
            }
            if (plot.types["combined"]){
                f.combined <- f.combined*f.capt
            }
        }
        ## Contour due to estimated bearings.
        if (plot.types["ang"] | plot.types["combined"] & fit$fit.types["ang"]){
            f.ang <- ang.density(fit, i, mask)
            if (plot.types["ang"]){
                show.contour(mask, f.x*f.ang, cols$ang)
            }
            if (plot.types["combined"]){
                f.combined <- f.combined*f.ang
            }
        }
        ## Contour due to estimated distances.
        if (plot.types["dist"] | plot.types["combined"] & fit$fit.types["dist"]){
            f.dist <- dist.density(fit, i, mask, dists)
            if (plot.types["dist"]){
                show.contour(mask, f.x*f.dist, cols$dist)
            }
            if (plot.types["combined"]){
                f.combined <- f.combined*f.dist
            }
        }
        ## Combined contour.
        if (plot.types["combined"]){
            show.contour(mask, f.combined, cols$combined)
        }
    }
    text(fit$traps, col = "red", labels = 1:nrow(fit$traps))
}

## Helper to get stuff in the right form for contour().
show.contour <- function(mask, dens, col = "black"){
    ## Divide densities by normalising constant before plotting.
    a <- attr(mask, "a")*10000
    ## Not conversion of area to square metres.
    dens <- dens/(a*sum(dens))
    unique.x <- sort(unique(mask[, 1]))
    unique.y <- sort(unique(mask[, 2]))
    z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
    n.mask <- nrow(mask)
    for (i in 1:n.mask){
        x <- mask[i, 1]
        y <- mask[i, 2]
        index.x <- which(x == unique.x)
        index.y <- which(y == unique.y)
        z[index.x, index.y] <- dens[i]
    }
    contour(unique.x, unique.y, z, add = TRUE, col = col)
}

## Calculating density due to estimated bearings.
ang.density <- function(fit, id, mask){
    capt <- fit$capt$bincapt[id, ]
    ang.capt <- fit$capt$ang[id, capt == 1]
    kappa <- getpar(fit, "kappa")
    mask.bearings <- bearings(fit$traps[capt == 1, , drop = FALSE], mask)
    mask.dens <- matrix(0, nrow = sum(capt), ncol = nrow(mask))
    for (i in 1:sum(capt)){
        mask.dens[i, ] <- dvm(ang.capt[i], mu = mask.bearings[i, ], kappa = kappa)
    }
    ## Plotting arrows.
    trappos <- fit$traps[which(capt == 1), , drop = FALSE]
    arrowlength <- 3
    sinb <- sin(ang.capt)*arrowlength
    cosb <- cos(ang.capt)*arrowlength
    arrows(trappos[, 1], trappos[, 2], trappos[, 1] + sinb, trappos[, 2] + cosb,
           length = 0.1)
    ## Returning densities.
    colProds(mask.dens)
}

## Calculating density due to estimated distances.
dist.density <- function(fit, id, mask, dists){
    capt <- fit$capt$bincapt[id, ]
    dists <- dists[capt == 1, ]
    dist.capt <- fit$capt$dist[id, capt == 1]
    alpha <- getpar(fit, "alpha")
    mask.dens <- matrix(0, nrow = sum(capt), ncol = nrow(mask))
    betas <- alpha/dists
    for (i in 1:sum(capt)){
        mask.dens[i, ] <- dgamma(dist.capt[i], shape = alpha, rate = betas[i, ])
    }
    ## Returning densities.
    colProds(mask.dens)
}


