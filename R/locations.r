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
#' \code{"dist"}, \code{"toa"}, \code{"combined"}, and \code{"all"},
#' where \code{"combined"} combines all information types together,
#' and \code{"all"} plots all possible contour types. When signal
#' strength information is used in the model fit, then selecting
#' \code{"capt"} here will use this information.
#' @param xlim A numeric vector of length 2, giving the x coordinate range.
#' @param ylim A numeric vector of length 2, giving the y coordinate range.
#' @param cols A list with named components corresponding to each
#' contour type (i.e., a subset of \code{"capt"}, \code{"ang"},
#' \code{"dist"}, \code{"toa"}, and \code{"combined"}). Each component
#' provides the colour the associated contour type (e.g., using a
#' character string such as \code{"red"}, or a call to the function
#' \link[grDevices]{rgb}).
#' @param plot.arrows Logical, if \code{TRUE} arrows indicating the
#' estimated bearing to the individual are plotted from detectors at
#' which detections were made.
#' @param plot.circles Logical, if \code{TRUE} circles indicating the
#' estimated distance to the individual are plotted around detectors
#' at which detections were made.
#' @param mask A matrix with two columns. Each row provides Cartesian
#' coordinates for the location of a mask point. The function
#' \link[admbsecr]{create.mask} will return a suitable object. The
#' mask used to fit the model \code{fit} will be used by default; this
#' argument is usually used when estimated location contours need to
#' be plotted to a higher resolution than this.
#' @param add Logical, if \code{TRUE} contours will be added to an
#' existing plot.
locations <- function(fit, id, infotypes = "combined",
                      xlim = range(mask[, 1]),
                      ylim = range(mask[, 2]),
                      cols = list(combined = "black", capt = "purple",
                          ang = "green", dist = "brown", toa = "blue"),
                      plot.arrows = any(c("ang", "all") %in% infotypes),
                      plot.circles = any(c("dist", "all") %in% infotypes),
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
    if (infotypes == "all"){
        infotypes <- c(fit$infotypes, "combined")
    }
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
            if (plot.arrows){
                show.arrows(fit, i)
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
            if (plot.circles){
                show.circles(fit, i)
            }
        }
        ## Combined contour.
        if (plot.types["combined"]){
            show.contour(mask, f.combined, cols$combined)
        }
    }
    points(fit$traps, col = "red", pch = 4, lwd = 2)
    if (length(id) == 1){
        points(fit$traps[capt == 1, ], col = "red", cex = 2, lwd = 2)
    }
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

## Plots arrows on traps where a detection was made, showing estimated bearing.
show.arrows <- function(fit, id){
    xlim <- par("usr")[c(1, 2)]
    ylim <- par("usr")[c(3, 4)]
    arrow.length <- 0.05*min(c(diff(range(xlim)), diff(range(ylim))))
    capt <- fit$capt$bincapt[id, ]
    ang.capt <- fit$capt$ang[id, capt == 1]
    trappos <- fit$traps[which(capt == 1), , drop = FALSE]
    sinb <- sin(ang.capt)*arrow.length
    cosb <- cos(ang.capt)*arrow.length
    arrows(trappos[, 1], trappos[, 2], trappos[, 1] + sinb, trappos[, 2] + cosb,
           length = 0.1, col = "red", lwd = 2)
}

## Plots circles around traps where a detection was made, showing estimated distance.
show.circles <- function(fit, id){
    capt <- fit$capt$bincapt[id, ]
    dist.capt <- fit$capt$dist[id, capt == 1]
    trappos <- fit$traps[which(capt == 1), , drop = FALSE]
    for (i in 1:nrow(trappos)){
        centre <- trappos[i, ]
        radius <- dist.capt[i]
        circles(centre, radius, col = "red", lwd = 2)
    }
}

circles <- function(centre, radius, ...){
    angs <- seq(0, 2*pi, length.out = 100)
    xs <- centre[1] + sin(angs)*radius
    ys <- centre[2] + cos(angs)*radius
    lines(xs, ys, ...)
}
