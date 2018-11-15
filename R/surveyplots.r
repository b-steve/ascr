#' Plotting mask and trap layout
#'
#' Plots the mask points and trap locations used in a model fitted
#' with the function \link{fit.ascr}.
#'
#' @param session The session for which the mask point and trap
#'     locations are to be plotted.
#' @param ... Further arguments to be passed to \link{plot}.
#' @inheritParams locations
#'
#' @examples
#' show.survey(example$fits$simple.hn)
#'
#' @export
show.survey <- function(fit, session = 1, ...){
    plot(get.mask(fit, session), pch = ".", cex = 3, asp = 1, ...)
    points(get.traps(fit, session), pch = 16, col = "red")
}

#' Plotting the detection probability surface
#'
#' Plots the detection probability surface, based on trap locations
#' and the estimated detection function from a model fitted using
#' \link{fit.ascr}.
#'
#' @inheritParams locations
#' @param session For multisession models, the session from which the
#'     plot should be generated.
#' @param surface Logical, if \code{TRUE} a 3D detection surface is
#'     plotted over the mask point locations, otherwise a contour plot
#'     is shown.
#' @param mask A matrix with two columns. Each row provides Cartesian
#'     coordinates for the location of a mask point. The function
#'     \link[ascr]{create.mask} will return a suitable object. The
#'     mask used to fit the model \code{fit} will be used by default;
#'     this argument is usually used when estimated location contours
#'     need to be plotted to a higher resolution than this.
#' @param col The colour of the plotted contours.
#' @param levels A numeric vector giving the values to be associated
#'     with the plotted contours. Alternatively, this can be the
#'     character string "esa", which results in a contour
#'     encapsulating an area equal to the estimated effective sampling
#'     area.
#' @param show.labels Logical, if \code{TRUE}, contours are labelled
#'     with the appropriate probability.
#' @param ... Arguments to be passed to \link{persp}.
#'
#' @examples
#' show.detsurf(example$fits$simple.hn)
#'
#' @export
show.detsurf <- function(fit, session = 1, surface = TRUE, mask = NULL, col = "black", levels = NULL, xlim = NULL, ylim = NULL,
                         show.labels = TRUE, trap.col = "red", add = FALSE, ...){
    match.esa <- FALSE
    if (!surface){
        if (is.character(levels)){
            if (levels == "esa"){
                match.esa <- TRUE
            } else {
                stop("If argument 'levels' is a character string, it must be \"esa\"")
            }
        }
    }
    if (is.null(mask)){
        mask <- get.mask(fit, session)
    }
    p.det <- p.dot(fit, session, points = mask)
    traps <- get.traps(fit, session)
    unique.x <- sort(unique(mask[, 1]))
    unique.y <- sort(unique(mask[, 2]))
    z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
    n.mask <- nrow(mask)
    n.traps <- nrow(traps)
    for (i in 1:n.mask){
        x <- mask[i, 1]
        y <- mask[i, 2]
        index.x <- which(x == unique.x)
        index.y <- which(y == unique.y)
        z[index.x, index.y] <- p.det[i]
    }
    ## Plotting surface.
    if (surface){
        perspmat <- persp(x = unique.x, y = unique.y, z = z, zlim = c(0, 1),
                          zlab = "", xlab = "", ylab = "", shade = 0.75,
                          col = "lightblue", theta = 30, phi = 30, border = NA,
                          box = TRUE, axes = FALSE, ticktype = "detailed", ...)
        ## Plotting z-axis.
        y.range <- diff(range(mask[, 2]))
        tick.start <- trans3d(x = rep(min(mask[, 1]), 2), y = rep(min(mask[, 2]), 2),
                              z = c(0, 1), pmat = perspmat)
        tick.end <- trans3d(x = rep(min(mask[, 1]), 2), y = rep(min(mask[, 2]), 2) - 0.025*y.range,
                            z = c(0, 1), pmat = perspmat)
        segments(x0 = tick.start$x, y0 = tick.start$y, x1 = tick.end$x, y1 = tick.end$y)
        label.pos <- trans3d(x = rep(min(mask[, 1]), 2), y = rep(min(mask[, 2]), 2) - 0.05*y.range,
                             z = c(0, 1), pmat = perspmat)
        text(x = label.pos$x, y = label.pos$y, labels = c(0, 1), srt = 15)
        title.pos <- label.pos <- trans3d(x = min(mask[, 1]), y = min(mask[, 2]) - 0.1*y.range,
                                          z = c(0.6), pmat = perspmat)
        text(x = title.pos$x, y = title.pos$y, labels = "Detection probability", srt = 105)
        for (i in 1:n.traps){
            ds <- distances(traps[i, , drop = FALSE], mask)
            nearest.mpoint <- which(ds == min(ds))[1]
            nearest.z <- p.det[nearest.mpoint]
            points(trans3d(x = traps[i, 1], y = traps[i, 2], z = nearest.z,
                           pmat = perspmat), pch = 16, col = trap.col)
        }
    } else {
        if (!add){
            plot(mask, type = "n", xlim = xlim, ylim = ylim, asp = 1,
                 xlab = "", ylab = "")
            points(traps, col = trap.col, pch = 4, lwd = 2)
        }
        if (match.esa){
            mask.area <- attr(mask, "area")
            esa <- coef(fit, "esa")
            n.inside <- round(esa/mask.area)
            levels <- sort(z, decreasing = TRUE)[n.inside]
        } else if (is.null(levels)){
            levels <- pretty(range(z, finite = TRUE), 10)
        }
        contour(x = unique.x, y = unique.y, z = z, levels = levels,
                col = col, drawlabels = show.labels, add = TRUE)
    }
}

#' Plotting an estimated density surface.
#'
#' Plots density surface estimated by a model fitted with the function
#' \link{fit.ascr}.
#'
#' @param zlim A numberic vector of length 2, giving the range of
#'     density contours.
#' @inheritParams locations
#'
#' @export
#'
#' @examples
#' simple.capt <- example$capt[1]
#' cov.df <- data.frame(x = example$mask[, 1]/1000, y = example$mask[, 2]/1000)
#' fit <- fit.ascr(capt = simple.capt, traps = example$traps, mask = example$mask,
#'                 fix = list(g0 = 1), ihd.opts = list(model = ~ x + y,
#'                                                     covariates = cov.df))
#' show.Dsurf(fit)
show.Dsurf <- function(fit, session = 1, xlim = NULL, ylim = NULL, zlim = NULL){
    D.mask <- fit$D.mask[[session]]
    mask <- get.mask(fit, session)
    if (is.null(xlim)){
        xlim <- range(mask[, 1])
    }
    if (is.null(ylim)){
        ylim <- range(mask[, 2])
    }
    mask.keep <- xlim[1] <= mask[, 1] & xlim[2] >= mask[, 1] &
        ylim[1] <= mask[, 2] & ylim[2] >= mask[, 2]
    mask <- mask[mask.keep, ]
    traps <- get.traps(fit, session)
    unique.x <- sort(unique(mask[, 1]))
    unique.y <- sort(unique(mask[, 2]))
    z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
    n.mask <- nrow(mask)
    n.traps <- nrow(traps)
    for (i in 1:n.mask){
        x <- mask[i, 1]
        y <- mask[i, 2]
        index.x <- which(x == unique.x)
        index.y <- which(y == unique.y)
        z[index.x, index.y] <- D.mask[mask.keep][i]
    }
    if (is.null(zlim)){
        zlim <- c(0, max(z, na.rm = TRUE))
    }
    levels <- pretty(zlim, 10)
    plot(mask, type = "n", asp = 1, xlab = "", ylab = "")
    image(x = unique.x, y = unique.y, z = z, zlim = zlim, add = TRUE)
    points(traps, col = "black", pch = 4, lwd = 2)
    contour(x = unique.x, y = unique.y, z = z, levels = levels,
            drawlabels = TRUE, add = TRUE)
}

#' Testing the mask object for a first calls model.
#'
#' Creates a diagnostic plot that can be used to test the adequacy of
#' a mask for a first calls model.
#'
#' @inheritParams fit.ascr
#' @inheritParams locations
#' @param pars A named list. Component names are parameter names, and
#'     components are values of parameters at which the mask is to be
#'     tested. Parameters must include \code{b0.ss}, \code{b1.ss}, and
#'     \code{sigma.ss}; see \link{fit.ascr} for further details.
#' @param cutoff The upper cutoff value; see \link{fit.ascr} for
#'     further details.
#' @param lower.cutoff The lower cutoff value; see \link{fit.ascr} for
#'     further details.
#' @param surface Logical, if \code{TRUE} a 3D detection surface is
#'     plotted over the mask point locations.
#' @param col Colours of points. Three colours are required if
#'     \code{surface} is \code{FALSE}; otherwise, just one is
#'     required.
#' @param pincol Colour of pins when \code{surface} is \code{TRUE}.
#' @param ... Arguments to be passed to \link{plot}.
#'
#' @export
mask.test <- function(fit = NULL, mask, traps, pars, cutoff, lower.cutoff, surface = FALSE,
                      col = c("black", "blue", "red"), pincol = "red", ...){
    if (!is.null(fit)){
        mask <- fit$args$mask
        traps <- fit$args$traps
        pars <- get.par(fit, pars = "fitted", as.list = TRUE)
        cutoff <- fit$args$ss.opts$cutoff
        lower.cutoff <- fit$args$ss.opts$lower.cutoff
    }
    traps.centroid <- matrix(apply(traps, 2, mean), nrow = 1)
    ss.means <- pars$b0.ss - pars$b1.ss*distances(mask, traps)
    log.p.au <- pnorm(cutoff, ss.means, pars$sigma.ss, lower.tail = FALSE, log.p = TRUE)
    log.p.al <- pnorm(lower.cutoff, ss.means, pars$sigma.ss, lower.tail = FALSE, log.p = TRUE)
    log.p.bu <- pnorm(cutoff, ss.means, pars$sigma.ss, lower.tail = TRUE, log.p = TRUE)
    log.p.bl <- pnorm(lower.cutoff, ss.means, pars$sigma.ss, lower.tail = TRUE, log.p = TRUE)
    n.mask <- nrow(mask)
    n.traps <- nrow(traps)
    n.combins <- 2^n.traps
    combins <- matrix(NA, nrow = n.combins, ncol = n.traps)
    for (i in 1:n.traps){
        combins[, i] <- rep(rep(c(0, 1), each = 2^(n.traps - i)), times = 2^(i - 1))
    }
    p <- numeric(n.mask)
    s <- numeric(n.mask)
    for (i in 1:n.mask){
        ## Log probabilities of detection and evasion at upper threshold.
        log.det.u.mat <- matrix(log.p.au[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
        log.evade.u.mat <- matrix(log.p.bu[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
        log.prob.u.mat <- log.det.u.mat
        log.prob.u.mat[combins == 0] <- log.evade.u.mat[combins == 0]
        ## Probabilities for each capture history.
        log.AU <- log(sum(exp(apply(log.prob.u.mat[-1, ], 1, sum))))
        log.BU <- sum(log.prob.u.mat[1, ])
        ## Log probabilities of detection and evasion at lower threshold.
        log.det.l.mat <- matrix(log.p.al[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
        log.evade.l.mat <- matrix(log.p.bl[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
        log.prob.l.mat <- log.det.l.mat
        log.prob.l.mat[combins == 0] <- log.evade.l.mat[combins == 0]
        ## Probabilities for each capture history.
        log.AL <- log(sum(exp(apply(log.prob.l.mat[-1, ], 1, sum))))
        log.BL <- sum(log.prob.l.mat[1, ])
        p[i] <- exp(log.AU)
        s[i] <- exp(log.AU + log.BL - log.AL)
    }
    if (surface){
        unique.x <- sort(unique(mask[, 1]))
        unique.y <- sort(unique(mask[, 2]))
        z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
        n.mask <- nrow(mask)
        for (i in 1:n.mask){
            x <- mask[i, 1]
            y <- mask[i, 2]
            index.x <- which(x == unique.x)
            index.y <- which(y == unique.y)
            z[index.x, index.y] <- s[i]
        }
        if (missing(col)){
            col = "lightblue"
        }
        perspmat <- persp(x = unique.x, y = unique.y, z = z, zlim = c(0, 1),
                          zlab = "", xlab = "", ylab = "", shade = 0.75,
                          col = col, theta = 30, phi = 30, border = NA,
                          box = TRUE, axes = FALSE, ticktype = "detailed", ...)
        ## Plotting z-axis.
        y.range <- diff(range(mask[, 2]))
        tick.start <- trans3d(x = rep(min(mask[, 1]), 2), y = rep(min(mask[, 2]), 2),
                              z = c(0, 1), pmat = perspmat)
        tick.end <- trans3d(x = rep(min(mask[, 1]), 2), y = rep(min(mask[, 2]), 2) - 0.025*y.range,
                            z = c(0, 1), pmat = perspmat)
        segments(x0 = tick.start$x, y0 = tick.start$y, x1 = tick.end$x, y1 = tick.end$y)
        label.pos <- trans3d(x = rep(min(mask[, 1]), 2), y = rep(min(mask[, 2]), 2) - 0.05*y.range,
                             z = c(0, 1), pmat = perspmat)
        text(x = label.pos$x, y = label.pos$y, labels = c(0, 1), srt = 15)
        title.pos <- label.pos <- trans3d(x = min(mask[, 1]), y = min(mask[, 2]) - 0.1*y.range,
                                          z = c(0.6), pmat = perspmat)
        text(x = title.pos$x, y = title.pos$y, labels = "", srt = 105)
        for (i in 1:n.traps){
            lines(trans3d(x = rep(traps[i, 1], 2), y = rep(traps[i, 2], 2), z = c(0, 0.1),
                          pmat = perspmat))
            points(trans3d(x = traps[i, 1], y = traps[i, 2], z = 0.1,
                          pmat = perspmat), pch = 16, col = pincol, cex = 0.5)
        }
    } else {
        plot(distances(mask, traps.centroid), s + p, type = "n", ylim = c(0, 1), ...)
        abline(h = c(0, 1), col = "lightgrey")
        points(distances(mask, traps.centroid), s, col = col[3], cex = 0.5)
        points(distances(mask, traps.centroid), p, col = col[2], cex = 0.5)
        points(distances(mask, traps.centroid), s + p, col = col[1], cex = 0.5)
    }
}


