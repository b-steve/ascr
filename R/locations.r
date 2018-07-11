#' Plotting estimated locations
#'
#' Plots estimated densities of animal locations, which are latent
#' variables in SECR models.
#'
#' @param fit A fitted model from \link{fit.ascr}.
#' @param id A numeric vector with row numbers from
#'     \code{fit$args$capt}, indicating which individuals' locations
#'     are to be plotted.
#' @param session The session with the detector array and invidual(s)
#'     to be plotted (for multi-session models only).
#' @param infotypes A character vector indicating the type(s) of
#'     information to be used when plotting the estimated density of
#'     location.  Elements can be a subset of \code{"capt"},
#'     \code{"bearing"}, \code{"dist"}, \code{"ss"}, \code{"toa"},
#'     \code{"combined"}, and \code{"all"}, where \code{"capt"} shows
#'     estimated location only using detection locations,
#'     \code{"combined"} combines all information types together, and
#'     \code{"all"} plots all possible contour types. When signal
#'     strength information is used in the model fit, \code{"capt"}
#'     and \code{"ss"} are equivalent as the signal strength
#'     information is built into the detection function. By default,
#'     only the most informative contour is plotted, i.e.,
#'     \code{"capt"} if the model was fitted with no additional
#'     information, and \code{"combined"} otherwise.
#' @param combine Logical, if \code{TRUE} then the information types
#'     specified in \code{infotypes} are combined into a single
#'     contour. If \code{FALSE} then separate contours are plotted for
#'     each information type.
#' @param xlim A numeric vector of length 2, giving the x coordinate
#'     range.
#' @param ylim A numeric vector of length 2, giving the y coordinate
#'     range.
#' @param mask A matrix with two columns. Each row provides Cartesian
#'     coordinates for the location of a mask point. The function
#'     \link[ascr]{create.mask} will return a suitable object. The
#'     mask used to fit the model \code{fit} will be used by default;
#'     this argument is usually used when estimated location contours
#'     need to be plotted to a higher resolution than this.
#' @param levels A numeric vector giving the values to be associated
#'     with the plotted contours.
#' @param nlevels The number of contour levels desired. Ignored if
#'     \code{levels} is provided.
#' @param density Logical, if \code{TRUE}, the labels on contours (and
#'     the levels specified by \code{levels}) refer to the density of
#'     the estimated distribution of the individual's location. If
#'     \code{FALSE}, the labels on contours (and the levels specified
#'     by \code{levels}) refer to the probability of the individual
#'     being located within the associated contour under the estimated
#'     distribution of the individual's location.
#' @param cols A list with named components corresponding to each
#'     contour type (i.e., a subset of \code{"capt"},
#'     \code{"bearing"}, \code{"dist"}, \code{"toa"}, and
#'     \code{"combined"}). Each component provides the colour of the
#'     associated contour type (e.g., using a character string such as
#'     \code{"red"}, or a call to the function
#'     \link[grDevices]{rgb}). By default, if only one contour is to
#'     be plotted, it will be plotted in black. Alternatively, a
#'     vector with a single element, specifying the colour for all
#'     contours.
#' @param ltys The line type of the contours, with the same required
#'     syntax as \code{cols}; see \link{par}.
#' @param trap.col The colour of the points representing detector
#'     locations.
#' @param circle.traps Logical, if \code{TRUE} circles are plotted
#'     around traps that made a detection of the individual in
#'     question.
#' @param show.labels Logical, if \code{TRUE}, contours are labelled
#'     with the appropriate probability density (if \code{density} is
#'     \code{TRUE}), or the corresponding probability of the
#'     individual being within the associated contour, under the
#'     estimated density (if \code{density} is \code{FALSE}).
#' @param plot.contours Logical, if \code{TRUE}, contours are
#'     plotted. Note that, if \code{FALSE}, nothing corresponding to
#'     the density of the individuals' locations is plotted unless
#'     \code{plot.estlocs} is \code{TRUE}.
#' @param plot.estlocs Logical, if \code{TRUE}, dots are plotted at
#'     the mode of the combined densities. If a density has more than
#'     a single mode (and the modes have the same density value) then
#'     a dot will be plotted for each.
#' @param keep.estlocs Logical, if \code{TRUE}, the locations of the
#'     estimated locations are returned.
#' @param plot.arrows Logical, if \code{TRUE}, arrows indicating the
#'     estimated bearing to the individual are plotted from detectors
#'     at which detections were made.
#' @param plot.circles Logical, if \code{TRUE}, circles indicating the
#'     estimated distance to the individual are plotted around
#'     detectors at which detections were made.
#' @param arrow.length Numeric, providing the length of the arrows
#'     (only used if \code{plot.arrows} is \code{TRUE}).
#' @param show.legend Logical, if \code{TRUE}, a legend will be added
#'     to the plot.
#' @param show.axes Logical, if \code{TRUE}, axes will be added to the
#'     plot.
#' @param add Logical, if \code{TRUE}, contours will be added to an
#'     existing plot.
#'
#' @examples
#' locations(example$fits$simple.hn, 1)
#' locations(example$fits$simple.hn, 1, levels = c(0.50, 0.90, 0.95))
#' \dontrun{
#' fine.mask <- create.mask(example$traps, 20, spacing = 0.2)
#' locations(example$fits$bearing.hn, 1, infotypes = "all", mask = fine.mask)
#' }
#'
#' @export
locations <- function(fit, id, session = 1, infotypes = NULL, combine = FALSE,
                      xlim = NULL, ylim = NULL, mask = get.mask(fit, session),
                      levels = NULL, nlevels = 10, density = FALSE,
                      cols = list(combined = "black", capt = "purple",
                          ss = "orange", bearing = "green", dist = "brown", toa = "blue"),
                      ltys = list(combined = "solid", capt = "solid",
                          ss = "solid", bearing = "solid", dist = "solid", toa = "solid"),
                      trap.col = "red", circle.traps = TRUE,
                      show.labels = TRUE, plot.contours = TRUE,
                      plot.estlocs = FALSE,
                      keep.estlocs = FALSE,
                      plot.arrows = "bearing" %in% fit$infotypes,
                      plot.circles = "dist" %in% fit$infotypes & !("bearing" %in% fit$infotypes),
                      arrow.length = NULL,
                      show.legend = FALSE, show.axes = TRUE, add = FALSE){
    ## Error if session argument is too large.
    if (session > fit$n.sessions){
        if (fit$n.sessions == 1){
            stop(paste("Argument 'session' too large; there was only ", fit$n.sessions, " session.", sep = ""))
        } else {
            stop(paste("Argument 'session' too large; there were only ", fit$n.sessions, " sessions.", sep = ""))
        }
    }
    ## Error for locations() with a directional model.
    if (!is.null(fit$args$ss.opts$directional)){
        if (fit$args$ss.opts$directional){
            stop("The locations() function has not yet been implemented for directional model fits.")
        }
    }
    if (!is.null(fit$args$ss.opts$het.source)){
        if (fit$args$ss.opts$het.source){
            stop("The locations() function has not yet been implemented for heterogeneous source strength model fits.")
        }
    }
    if (fit$first.calls){
        stop("The locations() function has not yet been implemented for first-call models.")
    }
    ## Error if combine specified without infotypes.
    if (missing(infotypes) & combine){
        stop("Argument `combine' is only useful if `infotypes' is provided.")
    }
    ## Saving estimated locations.
    if (keep.estlocs){
        estlocs <- matrix(0, nrow = length(id), ncol = 2)
        j <- 1
    }
    ## Sorting out limits.
    if (is.null(xlim)){
        xlim <- range(mask[, 1])
    }
    if (is.null(ylim)){
        ylim <- range(mask[, 2])
    }
    ## Setting up plotting area.
    if (!add){
        plot.new()
        plot.window(xlim = xlim, ylim = ylim, asp = 1)
        box()
        if (show.axes){
            axis(1)
            axis(2)
        }
    }
    ## Ignoring 'nlevels' if 'levels' is provided.
    if (!is.null(levels)){
        nlevels <- length(levels)
    }
    ## Logical value for whether or not any additional information was
    ## used in model fit.
    any.infotypes <- length(fit$infotypes) > 0
    ## Setting default infotypes.
    if (is.null(infotypes)){
        if (any.infotypes){
            infotypes <- "combined"
        } else {
            infotypes <- "capt"
        }
    }
    ## Error if "combined" is used when there is no additional information.
    if ("combined" %in% infotypes & !any.infotypes){
        stop("No additional information used in model 'fit', so a \"combined\" contour cannot be plotted.")
    }
    ## Working out which contours to plot.
    if ("all" %in% infotypes){
        infotypes <- c(fit$infotypes, "capt", "combined"[any.infotypes])
    }
    ## If "ss" is an infotype, set to "capt". OR NOT. <<- Wait what the hell is going on in this comment.
    ##infotypes[infotypes == "ss"] <- "capt"
    infotypes <- unique(infotypes)
    ## Setting colour to "black" if there is only one contour to be plotted.
    if (missing(cols)){
        if (length(infotypes) == 1){
            cols <- "black"
        }
    }
    if (missing(ltys)){
        if (length(infotypes) == 1){
            ltys <- "solid"
        }
    }
    if (length(cols) == 1){
        cols.save <- cols
        cols <- vector(mode = "list", length = length(infotypes))
        names(cols) <- infotypes
        cols[infotypes] <- cols.save
    }
    if (length(ltys) == 1){
        ltys.save <- ltys
        ltys <- vector(mode = "list", length = length(infotypes))
        names(ltys) <- infotypes
        ltys[infotypes] <- ltys.save
        if (combine){
            ltys["combined"] <- ltys.save
        }
    }
    plot.types <- c("combined", "capt", "ss", "bearing", "dist", "toa") %in% infotypes
    names(plot.types) <- c("combined", "capt", "ss", "bearing", "dist", "toa")
    if (combine){
        plot.types["combined"] <- TRUE
    }
    ## Setting all to TRUE if combined is used.
    ## Some error catching.
    for (i in c("bearing", "dist", "toa")){
        if (plot.types[i] & !fit$fit.types[i]){
            msg <- paste("Contours for information type '", i, "' cannot be plotted as this information was not used in the model 'fit'", sep = "")
            warning(msg)
            plot.types[i] <- FALSE
        }
    }
    traps <- get.traps(fit, session)
    detfn <- fit$args$detfn
    ss.link <- fit$args$ss.opts$ss.link
    dists <- distances(traps, mask)
    ## Calculating density due to animal locations.
    p.det <- p.dot(fit = fit, session = session, points = mask)
    ## Divide by normalising constant; not conversion to square metres.
    a <- attr(mask, "area")
    if (fit$fit.ihd){
        f.x <- p.det*fit$D.mask[[session]]/(a*sum(p.det*fit$D.mask[[session]]))
    } else {
        f.x <- p.det/(a*sum(p.det))
    }
    ## Extracting the session's capture history.
    capt.all <- get.capt(fit, session)
    ## Calculating conditional density of capture history, given location.
    for (i in id){
        if (plot.types["combined"]){
            if ((!combine) | (combine & plot.types["capt"])){
                f.combined <- f.x
            } else {
                f.combined <- 0*f.x + 1
            }
        }
        capt <- capt.all$bincapt[i, ]
        ## Contour due to capture history.
        if (plot.types["capt"] | plot.types["combined"] | plot.types["ss"]){
            det.pars <- get.par(fit, fit$detpars, as.list = TRUE)
            if (fit$fit.types["ss"]){
                det.pars$cutoff <- fit$args$ss.opts$cutoff
            }
            det.probs <- calc.detfn(dists, detfn, det.pars, ss.link)
            f.capt <- colProds(det.probs*capt + (1 - det.probs)*(1 - capt))
            if (plot.types["capt"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.capt, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$capt,
                                 lty = ltys$capt, show.labels = show.labels,
                                 plot.contours = plot.contours)
                }
            }
            if (fit$fit.types["ss"]){
                f.ss.capt <- ss.density(fit, i, session, mask, dists)
                f.ss <- f.ss.capt/f.capt
                ## Such a hack, but this keeps f.combined correct,
                ## below.
                f.capt <- f.ss.capt
                ## Ughhhh sorry about this one.
                f.ss[f.ss == Inf] <- 0
                if (plot.types["ss"]){
                    if (!combine){
                        show.contour(mask = mask, dens = f.x*f.ss, levels = levels,
                                     nlevels = nlevels, prob = !density, col = cols$ss,
                                     lty = ltys$ss, show.labels = show.labels,
                                     plot.contours = plot.contours)
                    }
                }
            } else {
                fit.ss <- f.capt*0 + 1
            }
            ## This shit is such a mess, sorry if you have to work out
            ## how this works later.
            if (plot.types["combined"] & !combine){
                f.combined <- f.combined*f.capt
            } else if (plot.types["combined"] & combine){
                f.true.capt <- f.capt/f.ss
                f.true.capt[f.ss == 0] <- 0
                if (plot.types["capt"] & plot.types["ss"]){
                    f.combined <- f.combined*f.capt
                } else if (plot.types["capt"] & !plot.types["ss"]){
                    f.combined <- f.combined*f.true.capt
                } else if (!plot.types["capt"] & plot.types["ss"]){
                    f.combined <- f.combined*f.ss
                }
            }
        }
        ## Contour due to estimated bearings.
        if (plot.types["bearing"] | plot.types["combined"] & fit$fit.types["bearing"]){
            f.bearing <- bearing.density(fit, i, session, mask)
            if (plot.types["bearing"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.bearing, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$bearing,
                                 lty = ltys$bearing, show.labels = show.labels,
                                 plot.contours = plot.contours)
                }
            }
            if (plot.types["combined"]){
                if ((!combine) | (combine & plot.types["bearing"])){
                    f.combined <- f.combined*f.bearing
                }
            }
        }
        ## Contour due to estimated distances.
        if (plot.types["dist"] | plot.types["combined"] & fit$fit.types["dist"]){
            f.dist <- dist.density(fit, i, session, mask, dists)
            if (plot.types["dist"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.dist, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$dist,
                                 lty = ltys$dist, show.labels = show.labels,
                                 plot.contours = plot.contours)
                }
            }
            if (plot.types["combined"]){
                if ((!combine) | (combine & plot.types["dist"])){
                    f.combined <- f.combined*f.dist
                }
            }
            if (plot.circles){
                show.circles(fit, i, session, trap.col)
            }
        }
        ## Contour due to measured times of arrival.
        if (plot.types["toa"] | plot.types["combined"] &
            fit$fit.types["toa"] & sum(capt) > 1){
            f.toa <- toa.density(fit, i, session, mask, dists)
            if (plot.types["toa"]){
                if (!combine){
                    show.contour(mask = mask, dens = f.x*f.toa, levels = levels,
                                 nlevels = nlevels, prob = !density, col = cols$toa,
                                 lty = ltys$toa, show.labels = show.labels,
                                 plot.contours = plot.contours)
                }
            }
            if (plot.types["combined"]){
                  if ((!combine) | (combine & plot.types["toa"])){
                      f.combined <- f.combined*f.toa
                  }
              }
        }
        ## Combined contour.
        if (plot.types["combined"]){
            show.contour(mask = mask, dens = f.combined, levels = levels,
                         nlevels = nlevels, prob = !density, col = cols$combined,
                         lty = ltys$combined, show.labels = show.labels,
                         plot.contours = plot.contours)
        }
        if (plot.estlocs){
            if ((any.infotypes & plot.types["combined"]) | !any.infotypes){
                f.estlocs <- if (any.infotypes) f.combined else f.capt*f.x
                mode.points <- which(f.estlocs == max(f.estlocs))
                points(mask[mode.points, 1], mask[mode.points, 2],
                       pch = 16, col = "black")
                if (keep.estlocs){
                    estlocs[j, ] <- c(mask[mode.points, 1], mask[mode.points, 2])
                    j <- j + 1
                }
            }
        }
    }
    ## Plotting traps, and circles around them.
    if (!add){
        points(traps, col = trap.col, pch = 4, lwd = 2)
        if (circle.traps){
            if (length(id) == 1){
                points(traps[capt == 1, , drop = FALSE], col = trap.col, cex = 2, lwd = 2)
            }
        }
    }
    ## Plotting arrows for estimated bearings.
    if (fit$fit.types["bearing"]){
        if (plot.arrows){
            if (fit$fit.types["dist"]){
                arrow.length <- capt.all$dist[i, capt == 1]
            }
            show.arrows(fit, i, session, arrow.length, trap.col)
        }
    }
    ## Making legend.
    if (show.legend){
        legend.labels <- infotypes
        legend.cols <- c(cols[infotypes], recursive = TRUE)
        legend.ltys <- c(ltys[infotypes], recursive = TRUE)
        legend("topright", legend = infotypes, lty = legend.ltys, col = legend.cols, bg = "white")
    }
    if (keep.estlocs){
        out <- list(estlocs = estlocs)
    } else {
        out <- invisible(TRUE)
    }
    out
}

## Helper to get stuff in the right form for contour().
show.contour <- function(mask, dens, nlevels, levels, prob, col = "black", lty = 1, show.labels, plot.contours){
    if (plot.contours){
        ## Divide densities by normalising constant before plotting.
        a <- attr(mask, "a")*10000
        ## Note conversion of area to square metres.
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
        ## Sorting out levels.
        if (is.null(levels)){
            levels <- pretty(range(z, finite = TRUE), nlevels)
        } else {
            if (prob){
                z.sort <- sort(z, decreasing = TRUE)
                probs.sort <- cumsum(z.sort)/sum(z.sort)
                prob.levels <- levels
                levels <- numeric(nlevels)
                for (i in 1:nlevels){
                    levels[i] <- z.sort[which(abs(probs.sort - prob.levels[i]) ==
                                              min(abs(probs.sort - prob.levels[i])))[1]]
                }
            }
        }
        if (prob){
            labels <- character(nlevels)
            for (i in 1:nlevels){
                labels[i] <- format(round(sum(z[z > levels[i]], na.rm = TRUE)/
                                          sum(z, na.rm = TRUE), 2), nsmall = 2)
            }
        } else {
            labels <- NULL
        }
        contour(x = unique.x, y = unique.y, z = z, levels = levels, labels = labels,
                col = col, lty = lty, drawlabels = show.labels, add = TRUE)
    }
}

## Calculating density due to estimated bearings.
bearing.density <- function(fit, id, session, mask){
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    bearing.capt <- capt.all$bearing[id, capt == 1]
    kappa <- get.par(fit, "kappa")
    mask.bearings <- bearings(get.traps(fit, session)[capt == 1, , drop = FALSE], mask)
    mask.dens <- matrix(0, nrow = sum(capt), ncol = nrow(mask))
    for (i in 1:sum(capt)){
        mask.dens[i, ] <- dvm(bearing.capt[i], mu = mask.bearings[i, ], kappa = kappa)
    }
    ## Returning densities.
    colProds(mask.dens)
}

## Calculating density due to estimated distances.
dist.density <- function(fit, id, session, mask, dists){
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    dists <- dists[capt == 1, , drop = FALSE]
    dist.capt <- capt.all$dist[id, capt == 1]
    alpha <- get.par(fit, "alpha")
    mask.dens <- matrix(0, nrow = sum(capt), ncol = nrow(mask))
    betas <- alpha/dists
    for (i in 1:sum(capt)){
        mask.dens[i, ] <- dgamma(dist.capt[i], shape = alpha, rate = betas[i, ])
    }
    ## Returning densities.
    colProds(mask.dens)
}

ss.density <- function(fit, id, session, mask, dists){
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    ss.capt <- capt.all$ss[id, ]
    det.pars <- get.par(fit, fit$detpars, cutoff = TRUE, as.list = TRUE)
    detfn <- fit$args$detfn
    ss.link <- fit$args$ss.opts$ss.link
    n.traps <- nrow(get.traps(fit, session))
    mask.dens <- matrix(0, nrow = n.traps, ncol = nrow(mask))
    for (i in 1:n.traps){
        if (capt[i] == 0){
            mask.dens[i, ] <- 1 - calc.detfn(dists[i, ], detfn, det.pars, ss.link)
        } else if (capt[i] == 1){
            mu.ss <- det.pars[["b0.ss"]] - det.pars[["b1.ss"]]*dists[i, ]
            mask.dens[i, ] <- dnorm(ss.capt[i], mu.ss, det.pars[["sigma.ss"]])
        } else {
            stop("The binary capture history must only contain 0s and 1s.")
        }
    }
    colProds(mask.dens)
}

toa.density <- function(fit, id, session, mask, dists){
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    dists <- dists[capt == 1, ]
    toa.capt <- capt.all$toa[id, capt == 1]
    sigma.toa <- get.par(fit, "sigma.toa")
    prod.times <- toa.capt - dists/fit$args$sound.speed
    toa.ssq <- aaply(prod.times, 2, function(x) sum((x - mean(x))^2))
    out <- (2*pi*sigma.toa^2)^((1 - sum(capt))/2)*
        exp(toa.ssq/(-2*sigma.toa^2))
}

## Plots arrows on traps where a detection was made, showing estimated bearing.
show.arrows <- function(fit, id, session, arrow.length = NULL, trap.col){
    xlim <- par("usr")[c(1, 2)]
    ylim <- par("usr")[c(3, 4)]
    if (is.null(arrow.length)){
        arrow.length <- 0.05*min(c(diff(range(xlim)), diff(range(ylim))))
    }
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    bearing.capt <- capt.all$bearing[id, capt == 1]
    trappos <- get.traps(fit, session)[which(capt == 1), , drop = FALSE]
    sinb <- sin(bearing.capt)*arrow.length
    cosb <- cos(bearing.capt)*arrow.length
    arrows(trappos[, 1], trappos[, 2], trappos[, 1] + sinb, trappos[, 2] + cosb,
           length = 0.1, col = trap.col, lwd = 2)
}

## Plots circles around traps where a detection was made, showing estimated distance.
show.circles <- function(fit, id, session, trap.col){
    capt.all <- get.capt(fit, session)
    capt <- capt.all$bincapt[id, ]
    dist.capt <- capt.all$dist[id, capt == 1]
    trappos <- get.traps(fit, session)[which(capt == 1), , drop = FALSE]
    for (i in 1:nrow(trappos)){
        centre <- trappos[i, ]
        radius <- dist.capt[i]
        circles(centre, radius, col = trap.col, lwd = 2)
    }
}

circles <- function(centre, radius, ...){
    bearings <- seq(0, 2*pi, length.out = 100)
    xs <- centre[1] + sin(bearings)*radius
    ys <- centre[2] + cos(bearings)*radius
    lines(xs, ys, ...)
}
