#' Plotting estimated animal location densities.
#'
#' Plots a density of individual animals' locations from a fit returned by
#' \code{admbsecr()}.
#'
#' @param fit a fitted model returned by \code{\link[admbsecr]{admbsecr}}.
#' @param ... additional arguments to be passed to
#' \code{\link[graphics]{contour}}.
#' @rdname contours
#' @export
contours <- function(fit, ...){
  UseMethod("contours", fit)
}

#' @S3method contours default
contours.default <- function(fit, ...){
  stop(paste("This function does not work with admbsecr fits of method ",
             "\"", fit$method, "\"", sep = ""))
}

#' @rdname contours
#' @param dets which individuals' location densities are plotted.
#' @param add logical, if \code{TRUE} the contours are added to an already
#' existing plot.
#' @param heat logical, if \code{TRUE} a levelplot is used instead of contours.
#' @param cols a vector specifying the colours of the main, extra, and
#' simple contours, in that order. Can be of length 1 if all three are
#' to be the same.
#' @param ltys vector specifying the line types for the main, simple,
#' and extra contours, in that order. Can be of length 1 if all three
#' are to be the same.
#' @param trapnos logical, if \code{TRUE} the trap identification numbers are
#' plotted.
#' @param problevels a vector indicating which probabilities should be
#' associated with the levels of the contours.
#' @param showcapt logical, if \code{TRUE} circles are drawn around detectors
#' on which the detection was made.
#' @param xlim numeric vector of length 2, giving the x coordinate range.
#' @param ylim numeric vector of length 2, giving the y coordinate range.
#' @method contours simple
#' @S3method contours simple
contours.simple <- function(fit, dets = "all", add = FALSE, heat = FALSE,
                            cols = "black", ltys = 1, trapnos = FALSE,
                            problevels = NULL,
                            showcapt = length(dets) == 1 && dets != "all" && !add,
                            xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols,
                                        ltys, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  mask <- fit$mask
  allcapt <- data$capt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  if (fit$detfn == "hn"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    allprobs <- g0*exp(-dist^2/(2*sigma^2))
  } else if (fit$detfn == "hr"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    z <- parvals["z"]
    allprobs <- g0*(1 - exp(-(dist/sigma)^(-z)))
  } else if (fit$detfn == "th"){
    shape <- parvals["shape"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape - scale*dist)
  } else if (fit$detfn == "logth"){
    shape1 <- parvals["shape1"]
    shape2 <- parvals["shape2"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape1 - exp(shape2 + scale*dist))
  }
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    maskdens <- exp(simpledens)*D
    maskdens <- maskdens/sum(maskdens)
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}

#' @rdname contours
#' @param partition if logical, \code{TRUE} indicates that the contributions
#' to the countour due to both the binary capture history data and the
#' supplementary information are also to be plotted. If character, \code{"all"}
#' and \code{"none"} correspond to \code{"TRUE"} and \code{"FALSE"}
#' respectively. \code{"simple"} indicates only to add contours due to the
#' binary capture history data, and \code{"extra"} indicates only to add
#' contours due to the supplementary information.
#' @method contours toa
#' @S3method contours toa
contours.toa <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                         heat = FALSE, cols = c("black", rgb(0, 1, 0, 0.4),
                                         rgb(0, 0, 1, 0.4)), ltys = 1,
                         trapnos = FALSE, problevels = NULL,
                         showcapt = length(dets) == 1 && dets != "all" && !add,
                         xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols, ltys,
                                        partition = partition, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  extra.contours <- check.partition(partition)
  plot.simple <- extra.contours$plot.simple
  plot.extra <- extra.contours$plot.extra
  plot.part <- plot.simple | plot.extra
  mask <- fit$mask
  allcapt <- data$capt
  alltoacapt <- data$toacapt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  sigmatoa <- parvals["sigmatoa"]
  if (fit$detfn == "hn"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    allprobs <- g0*exp(-dist^2/(2*sigma^2))
  } else if (fit$detfn == "hr"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    z <- parvals["z"]
    allprobs <- g0*(1 - exp(-(dist/sigma)^(-z)))
  } else if (fit$detfn == "th"){
    shape <- parvals["shape"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape - scale*dist)
  } else if (fit$detfn == "logth"){
    shape1 <- parvals["shape1"]
    shape2 <- parvals["shape2"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape1 - exp(shape2 + scale*dist))
  }
  times <- dist/330
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    ## Can only incorporate TOA part if more than one detection.
    if (sum(allcapt[i, ]) > 1){
      toadens <- logdens.toa(alltoacapt, allcapt, times, sigmatoa, i)
    } else {
      if (partition){
        warning("Setting partition to FALSE; no TOA information.")
        plot.part <- FALSE
      }
      toadens <- 0
    }
    maskdens <- exp(simpledens + toadens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(simpledens, toadens, plot.simple, plot.extra, D, mask,
                          problevels, cols = cols[2:3], ltys = ltys[2:3], ...)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}

## No partition on signal strength as it is included in the detection function.
#' @rdname contours
#' @method contours ss
#' @S3method contours ss
contours.ss <- function(fit, dets = "all", add = FALSE, heat = FALSE,
                        cols = c("black", rgb(0, 1, 0, 0.4),
                                         rgb(0, 0, 1, 0.4)), ltys = 1,
                        trapnos = FALSE, problevels = NULL,
                        showcapt = length(dets) == 1 && dets != "all" && !add,
                        xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols,
                                        ltys, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  mask <- fit$mask
  allcapt <- data$capt
  allsscapt <- data$sscapt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  nmask <- data$nmask
  cutoff <- fit$data$c
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  ssb0 <- parvals["ssb0"]
  ssb1 <- parvals["ssb1"]
  sigmass <- parvals["sigmass"]
  lpred <- ssb0 + ssb1*dist
  if (fit$detfn == "identity") muss <- lpred else muss <- exp(lpred)
  allnonprobs <- pnorm(cutoff, muss, sigmass)
  for (i in dets){
    ssdens <- logdens.ss(allcapt, allsscapt, allnonprobs, ntraps,
                         muss, sigmass, i)
    maskdens <- exp(ssdens)*D
    maskdens <- maskdens/sum(maskdens)
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}

#' @rdname contours
#' @method contours sstoa
#' @S3method contours sstoa
contours.sstoa <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                           heat = FALSE, cols = c("black", rgb(0, 1, 0, 0.4),
                                           rgb(0, 0, 1, 0.4)), ltys = 1,
                           trapnos = FALSE, problevels = NULL,
                           showcapt = length(dets) == 1 && dets != "all" && !add,
                           xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols,
                                        ltys, partition = partition, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  extra.contours <- check.partition(partition)
  plot.simple <- extra.contours$plot.simple
  plot.extra <- extra.contours$plot.extra
  plot.part <- plot.simple | plot.extra
  mask <- fit$mask
  allcapt <- data$capt
  allsscapt <- data$sscapt
  alltoacapt <- data$toacapt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  nmask <- data$nmask
  cutoff <- fit$data$c
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  ssb0 <- parvals["ssb0"]
  ssb1 <- parvals["ssb1"]
  sigmass <- parvals["sigmass"]
  sigmatoa <- parvals["sigmatoa"]
  times <- dist/330
  muss <- ssb0 + ssb1*dist
  allnonprobs <- pnorm(cutoff, muss, sigmass)
  for (i in dets){
    ssdens <- logdens.ss(allcapt, allsscapt, allnonprobs, ntraps,
                         muss, sigmass, i)
    ## Can only incorporate TOA part if more than one detection.
    if (sum(allcapt[i, ]) > 1){
      toadens <- logdens.toa(alltoacapt, allcapt, times, sigmatoa, i)
    } else {
      if (partition){
        warning("Setting partition to FALSE; no TOA information.")
        plot.part <- FALSE
      }
      toadens <- 0
    }
    maskdens <- exp(ssdens + toadens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(ssdens, toadens, plot.simple, plot.extra, D, mask,
                          problevels, cols = cols[2:3], ltys = ltys[2:3], ...)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}


#' @rdname contours
#' @param arrows logical, indicating whether or not arrows are to be
#' plotted to show estimated animal bearing from traps.
#' @method contours ang
#' @S3method contours ang
contours.ang <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                         heat = FALSE, cols = c("black", rgb(0, 1, 0, 0.4),
                                         rgb(0, 0, 1, 0.4)), ltys = 1,
                         trapnos = FALSE, problevels = NULL,
                         showcapt = length(dets) == 1 && dets != "all" && !add,
                         arrows = showcapt, xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols, ltys,
                                        partition = partition, arrows = arrows, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
  arrows <- updated.arguments$arrows
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  extra.contours <- check.partition(partition)
  plot.simple <- extra.contours$plot.simple
  plot.extra <- extra.contours$plot.extra
  plot.part <- plot.simple | plot.extra
  mask <- fit$mask
  allcapt <- data$capt
  allangcapt <- data$angcapt
  traps <- fit$traps
  dist <- data$dist
  ang <- data$ang
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  kappa <- parvals["kappa"]
  if (fit$detfn == "hn"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    allprobs <- g0*exp(-dist^2/(2*sigma^2))
  } else if (fit$detfn == "hr"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    z <- parvals["z"]
    allprobs <- g0*(1 - exp(-(dist/sigma)^(-z)))
  } else if (fit$detfn == "th"){
    shape <- parvals["shape"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape - scale*dist)
  } else if (fit$detfn == "logth"){
    shape1 <- parvals["shape1"]
    shape2 <- parvals["shape2"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape1 - exp(shape2 + scale*dist))
  }
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    angdens <- logdens.ang(allangcapt, allcapt, ang, kappa, i)
    maskdens <- exp(simpledens + angdens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(simpledens, angdens, plot.simple, plot.extra, D, mask,
                          problevels, cols = cols[2:3], ltys = ltys[2:3], ...)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
  if (arrows){
    plot.arrows(traps, allcapt, allangcapt, i, heat)
  }
}

#' @rdname contours
#' @method contours disttc
#' @S3method contours disttc
contours.disttc <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                            heat = FALSE, cols = c("black", rgb(0, 1, 0, 0.4)),
                            ltys = 1, trapnos = FALSE, problevels = NULL,
                            showcapt = length(dets) == 1 && dets != "all" && !add,
                            xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols, ltys,
                                        partition = partition, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  extra.contours <- check.partition(partition)
  plot.simple <- extra.contours$plot.simple
  plot.extra <- extra.contours$plot.extra
  plot.part <- plot.simple | plot.extra
  mask <- fit$mask
  allcapt <- data$capt
  alldistcapt <- data$distcapt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  alpha <- parvals["alpha"]
  g01 <- parvals["g01"]
  g02 <- parvals["g02"]
  sigma1 <- parvals["sigma1"]
  sigma2 <- parvals["sigma2"]
  allprobs <- matrix(0, nrow = nrow(dist), ncol = ncol(dist))
  allprobs[1, ] <- g01*exp(-dist[1, ]^2/(2*sigma1^2))
  allprobs[2, ] <- g02*exp(-dist[2, ]^2/(2*sigma2^2))
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    distdens <- logdens.dist(alldistcapt, allcapt, dist, alpha, i)
    maskdens <- exp(simpledens + distdens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(simpledens, distdens, plot.simple, plot.extra, D, mask,
                          problevels, cols = cols[2:3], ltys = ltys[2:3], ...)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
  plot.circles(traps, allcapt, alldistcapt, i, heat)
}

#' @rdname contours
#' @method contours dist
#' @S3method contours dist
contours.dist <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                          heat = FALSE, cols = c("black", rgb(0, 1, 0, 0.4)),
                          ltys = 1, trapnos = FALSE, problevels = NULL,
                          showcapt = length(dets) == 1 && dets != "all" && !add,
                          xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, cols, ltys,
                                        partition = partition, ...)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
  cols <- updated.arguments$cols
  ltys <- updated.arguments$ltys
  extra.contours <- check.partition(partition)
  plot.simple <- extra.contours$plot.simple
  plot.extra <- extra.contours$plot.extra
  plot.part <- plot.simple | plot.extra
  mask <- fit$mask
  allcapt <- data$capt
  alldistcapt <- data$distcapt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  allpars <- fit$parnames
  estpars <- names(coefs)
  parvals <- numeric(length(allpars))
  names(parvals) <- allpars
  for (i in allpars){
    parvals[i] <- ifelse(i %in% estpars, coefs[i], data[[i]])
  }
  D <- parvals["D"]
  alpha <- parvals["alpha"]
  if (fit$detfn == "hn"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    allprobs <- g0*exp(-dist^2/(2*sigma^2))
  } else if (fit$detfn == "hr"){
    g0 <- parvals["g0"]
    sigma <- parvals["sigma"]
    z <- parvals["z"]
    allprobs <- g0*(1 - exp(-(dist/sigma)^(-z)))
  } else if (fit$detfn == "th"){
    shape <- parvals["shape"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape - scale*dist)
  } else if (fit$detfn == "logth"){
    shape1 <- parvals["shape1"]
    shape2 <- parvals["shape2"]
    scale <- parvals["scale"]
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    allprobs <- 0.5 - 0.5*erf(shape1 - exp(shape2 + scale*dist))
  }
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    distdens <- logdens.dist(alldistcapt, allcapt, dist, alpha, i)
    maskdens <- exp(simpledens + distdens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(simpledens, distdens, plot.simple, plot.extra, D, mask,
                          problevels, cols = cols[2:3], ltys = ltys[2:3], ...)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat,
                      problevels, cols[1], ltys[1], ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
  plot.circles(traps, allcapt, alldistcapt, i, heat)
}

#' @rdname contours
#' @method contours mrds
#' @S3method contours mrds
contours.mrds <- function(fit, dets = "all", add = FALSE, trapnos = FALSE,
                          showcapt = length(dets) == 1 && dets != "all" && !add,
                          xlim = NULL, ylim = NULL, ...){
  if (length(dets) > 1 & showcapt){
    warning("Setting showcapt to FALSE as length(dets) > 1")
    showcapt <- FALSE
  }
  data <- fit$data
  n <- data$n
  mask <- fit$mask
  allcapt <- data$capt
  alldistcapt <- data$indivdist
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add){
    make.plot(mask, xlim, ylim)
  }
  for (i in dets){
    plot.traps(traps, allcapt, i, FALSE, trapnos, showcapt)
    plot.circles(traps, allcapt, alldistcapt, i, FALSE)
  }
}

## Checks inputs and returns altered argument values.
warning.contours <- function(n, dets, add, heat, showcapt, cols, ltys,
                             partition = NULL, arrows = NULL, ...){
  if (length(dets) == 1 && dets == "all"){
    dets <- 1:n
  }
  elist <- list(...)
  refusepars <- c("x", "y", "z", "levels", "labels", "col", "lty")
  badpars <- names(elist)[names(elist) %in% refusepars]
  if (length(badpars) == 1){
    stop(paste(badpars, "cannot be given as an argument to contour() via \"...\"."))
  } else if (length(badpars) > 1){
    stop(paste(paste(badpars, collapse = ", "), "cannot be given as arguments to contour() via \"...\"."))
  }
  if (add & heat){
    warning("Setting add to FALSE as heat is TRUE")
  }
  if (length(dets) > 1 & heat){
    stop("Only one animal can be plotted when heat is TRUE")
  }
  if (length(dets) > 1 & showcapt){
    warning("Setting showcapt to FALSE as length(dets) > 1")
    showcapt <- FALSE
  }
  if (add & showcapt){
    warning("Setting showcapt to FALSE as add is TRUE")
    showcapt <- FALSE
  }
  if (!is.null(partition)){
    if (heat & !(partition == FALSE | partition == "none")){
      warning("Setting partition to FALSE as heat is TRUE")
      partition <- FALSE
    }
    if (heat & !(partition == FALSE | partition == "none")){
      warning("Setting partition to FALSE as heat is TRUE")
      partition <- FALSE
    }
    if (length(dets) > 1 & !(partition == FALSE | partition == "none")){
      warning("Setting partition to FALSE as length(dets) > 1")
      partition <- FALSE
    }
  }
  if (length(cols) == 1){
    cols <- rep(cols, 3)
  }
  if (length(ltys) == 1){
    ltys <- rep(ltys, 3)
  }
  list(dets = dets, showcapt = showcapt, partition = partition, arrows = arrows,
       cols = cols, ltys = ltys)
}

## Sets up indicators for simple and extra contour plotting.
check.partition <- function(partition, cols, ltys){
  if (is.logical(partition)){
    plot.simple <- partition
    plot.extra <- partition
  } else {
    if (partition == "all"){
      plot.simple <- TRUE
      plot.extra <- TRUE
    } else if (partition == "none"){
      plot.simple <- FALSE
      plot.extra <- FALSE
    } else if (partition == "simple"){
      plot.simple <- TRUE
      plot.extra <- FALSE
    } else if (partition == "extra"){
      plot.simple <- FALSE
      plot.extra <- TRUE
    } else {
      stop("partition must be \"all\", \"none\", \"simple\" or \"extra\"")
    }
  }
  list(plot.simple = plot.simple, plot.extra = plot.extra)
}

## Generates the axes and plotting area.
make.plot <- function(mask, xlim, ylim){
  if (is.null(xlim)) xlim <- range(mask[, 1])
  if (is.null(ylim)) ylim <- range(mask[, 2])
  plot(mask, type = "n", xlim = xlim, ylim = ylim, asp = 1)
}

## Calculates the log of the animal density due to binary capture history data.
logdens.simple <- function(allcapt, allprobs, ntraps, i){
  capt <- allcapt[i, ]
  probs <- allprobs
  for (j in 1:ntraps){
    if (capt[j] == 0) probs[j, ] <- 1 - probs[j, ]
  }
  apply(log(probs), 2, sum)
}

## Calculates the log of the animal density due to signal strength data.
logdens.ss <- function(allcapt, allsscapt, allnonprobs, ntraps, muss, sigmass, i){
  capt <- allcapt[i, ]
  sscapt <- allsscapt[i, ]
  probs <- matrix(0, nrow = ntraps, ncol = ncol(allnonprobs))
  for (j in 1:ntraps){
    if (capt[j] == 1){
      probs[j, ] <- dnorm(sscapt[j], muss[j, ], sigmass)
    } else {
      probs[j, ] <- allnonprobs[j, ]
    }
  }
  apply(log(probs), 2, sum)
}

## Calculates the log of the animal density due to angle data.
logdens.ang <- function(allangcapt, allcapt, ang, kappa, i){
  capt <- allcapt[i, ]
  angcapt <- allangcapt[i, ]
  dettraps <- which(capt == 1)
  dens.ang <- dvm(angcapt[dettraps], ang[dettraps, ], kappa)
  if (length(dettraps) == 1){
    dim(dens.ang) <- c(1, ncol(ang))
  }
  apply(log(dens.ang), 2, sum)
}

## Calculates the log of the animal density due to TOA data.
logdens.toa <- function(alltoacapt, allcapt, times, sigmatoa, i){
  capt <- allcapt[i, ]
  toacapt <- alltoacapt[i, ]
  dettraps <- which(capt == 1)
  dettimes <- times[dettraps, ]
  toacapt <- toacapt[dettraps]
  esttimes <- toacapt - dettimes
  ssqtoa <- apply(esttimes, 2, function(x) sum((x - mean(x))^2))
  (1 - sum(capt))*log(sigmatoa^2) - ssqtoa/(2*sigmatoa^2)
}

logdens.dist <- function(alldistcapt, allcapt, dist, alpha, i){
  capt <- allcapt[i, ]
  distcapt <- alldistcapt[i, ]
  dettraps <- which(capt == 1)
  beta <- alpha/dist[dettraps, ]
  dens.dist <- dgamma(distcapt[dettraps], alpha, beta)
  if (length(dettraps) == 1){
    dim(dens.dist) <- c(1, ncol(dist))
  }
  apply(log(dens.dist), 2, sum)
}

## Plots the overall contour for the animal.
plot.main.contour <- function(maskdens, mask, xlim, ylim, heat,
                              problevels, col, lty, ...){
  x <- mask[, 1]
  y <- mask[, 2]
  unique.x <- sort(unique(x))
  unique.y <- sort(unique(y))
  z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
  for (j in 1:length(maskdens)){
    xind <- which(unique.x == x[j])
    yind <- which(unique.y == y[j])
    z[xind, yind] <- maskdens[j]
  }
  if (heat){
    if (is.null(xlim)) xlim <- range(x)
    if (is.null(ylim)) ylim <- range(y)
    image(x = unique.x, y = unique.y, z = z, xlab = "x", ylab = "y",
          xlim = xlim, ylim = ylim, asp = 1)
    box()
  } else {
    if (is.null(problevels)){
      nlevels <- list(...)$nlevels
      if (is.null(nlevels)) nlevels <- 10
      zlim <- range(z, na.rm = TRUE)
      levels <- seq(zlim[1], zlim[2], length.out = nlevels)
    } else {
      nlevels <- length(problevels)
      zs <- sort(z, decreasing = TRUE)
      probs <- cumsum(zs)/sum(zs)
      levels <- numeric(nlevels)
      for (i in 1:nlevels){
        levels[i] <- zs[which(abs(probs - problevels[i]) ==
                              min(abs(probs - problevels[i])))[1]]
      }
    }
    labels <- numeric(nlevels)
    for (i in 1:nlevels){
      labels[i] <- format(round(sum(z[z > levels[i]], na.rm = TRUE)/
                                sum(z, na.rm = TRUE), 2), nsmall = 2)
    }
    contour(x = unique.x, y = unique.y, z = z, add = TRUE, levels = levels,
            labels = labels, col = col, lty = lty, ...)
  }
}

## Plots the extra contours (i.e., when partition != FALSE).
plot.other.contours <- function(detdens, otherdens, plot.simple,
                                plot.extra, D, mask, problevels,
                                cols, ltys, ...){
  secr.maskdens <- exp(detdens)*D
  secr.maskdens <- secr.maskdens/sum(secr.maskdens)
  other.maskdens <- exp(otherdens)*D
  other.maskdens <- other.maskdens/sum(other.maskdens)
  x <- mask[, 1]
  y <- mask[, 2]
  unique.x <- sort(unique(x))
  unique.y <- sort(unique(y))
  z1 <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
  z2 <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
  for (j in 1:length(detdens)){
      xind <- which(unique.x == x[j])
      yind <- which(unique.y == y[j])
      z1[xind, yind] <- secr.maskdens[j]
      z2[xind, yind] <- other.maskdens[j]
  }
  if (plot.simple){
    if (is.null(problevels)){
      nlevels <- list(...)$nlevels
      if (is.null(nlevels)) nlevels <- 10
      zlim <- range(z1, na.rm = TRUE)
      levels <- seq(zlim[1], zlim[2], length.out = nlevels)
    } else {
      nlevels <- length(problevels)
      zs <- sort(z1, decreasing = TRUE)
      probs <- cumsum(zs)/sum(zs)
      levels <- numeric(nlevels)
      for (i in 1:nlevels){
        levels[i] <- zs[which(abs(probs - problevels[i]) ==
                              min(abs(probs - problevels[i])))[1]]
      }
    }
    labels <- numeric(nlevels)
    for (i in 1:nlevels){
      labels[i] <- format(round(sum(z1[z1 > levels[i]], na.rm = TRUE)/
                                sum(z1, na.rm = TRUE), 2), nsmall = 2)
    }
    contour(x = unique.x, y = unique.y, z = z1, add = TRUE, levels = levels,
            labels = labels, col = cols[2], lty = ltys[2], ...)
  }
  if (plot.extra){
    if (is.null(problevels)){
      nlevels <- list(...)$nlevels
      if (is.null(nlevels)) nlevels <- 10
      zlim <- range(z2, na.rm = TRUE)
      levels <- seq(zlim[1], zlim[2], length.out = nlevels)
    } else {
      nlevels <- length(problevels)
      zs <- sort(z2, decreasing = TRUE)
      probs <- cumsum(zs)/sum(zs)
      levels <- numeric(nlevels)
      for (i in 1:nlevels){
        levels[i] <- zs[which(abs(probs - problevels[i]) ==
                              min(abs(probs - problevels[i])))[1]]
      }
    }
    labels <- numeric(nlevels)
    for (i in 1:nlevels){
      labels[i] <- format(round(sum(z2[z2 > levels[i]], na.rm = TRUE)/
                                sum(z2, na.rm = TRUE), 2), nsmall = 2)
    }
    contour(x = unique.x, y = unique.y, z = z2, add = TRUE, levels = levels,
            labels = labels, col = cols[1], lty = ltys[1], ...)
  }
}

## Plots the traps.
plot.traps <- function(traps, allcapt, i, heat, trapnos, showcapt){
  trapcol <- ifelse(heat, "black", "red")
  if (trapnos){
    text(traps, labels = 1:nrow(traps), col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol)
  }
  if (showcapt){
    points(traps[which(allcapt[i, ] == 1), , drop = FALSE], cex = 2,
           lwd = 2, col = trapcol)
  }
}

## Plots arrows on traps for the angles method.
plot.arrows <- function(traps, allcapt, allangcapt, i, heat){
  xlim <- par("usr")[1:2]
  ylim <- par("usr")[3:4]
  arrowlength <- 0.15*min(range(xlim), range(ylim))
  arrowcol <- ifelse(heat, "black", "red")
  trappos <- traps[which(allcapt[i, ] == 1), , drop = FALSE]
  bearings <- allangcapt[i, which(allcapt[i, ] == 1)]
  sinb <- -sin(bearings)*arrowlength
  cosb <- -cos(bearings)*arrowlength
  arrows(trappos[, 1], trappos[, 2], trappos[, 1] + sinb, trappos[, 2] + cosb,
         length = 0.1, col = arrowcol)
}

## Plots circles indicating estimated distance from traps.
plot.circles <- function(traps, allcapt, alldistcapt, i, heat){
  circlecol <- ifelse(heat, "black", "red")
  trappos <- as.matrix(traps[which(allcapt[i, ] == 1), , drop = FALSE])
  dists <- alldistcapt[i, which(allcapt[i, ] == 1)]
  for (j in 1:nrow(trappos)){
    circles(trappos[j, ], dists[j], col = circlecol, lwd = 2)
  }
}

circles <- function(cent, rad, col = "black", lwd = 1){
  angs <- seq(0, 2*pi, length.out = 100)
  xs <- cent[1] + sin(angs)*rad
  ys <- cent[2] + cos(angs)*rad
  lines(xs, ys, col = col, lwd = lwd)
}
