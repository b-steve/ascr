#' Plotting estimated animal location densities.
#'
#' Plots a density of individual animals' locations from a fit returned by
#' \code{admbsecr()}.
#'
#' @param fit a fitted model returned by \code{\link[admbsecr]{admbsecr}}.
#' @param ... additional arguments to be passed to \code{\link[graphics]{contour}}.
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
#' @param col specifies the colour of the contours to be plotted.
#' @param trapnos logical, if \code{TRUE} the trap identification numbers are
#' plotted.
#' @param showcapt logical, if \code{TRUE} circles are drawn around detectors
#' on which the detection was made.
#' @param xlim numeric vector of length 2, giving the x coordinate range.
#' @param ylim numeric vector of length 2, giving the y coordinate range.
#' @method contours simple
#' @S3method contours simple
contours.simple <- function(fit, dets = "all", add = FALSE, heat = FALSE,
                            col = "black", trapnos = FALSE,
                            showcapt = length(dets) == 1 && dets != "all",
                            xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  mask <- fit$mask
  allcapt <- data$capt
  traps <- fit$traps
  dist <- data$dist
  ntraps <- data$ntraps
  coefs <- coef(fit)
  if (!add & !heat){
    make.plot(mask, xlim, ylim)
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    maskdens <- exp(simpledens)*D
    maskdens <- maskdens/sum(maskdens)
    plot.main.contour(maskdens, mask, xlim, ylim, heat, col, ...)
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
                         heat = FALSE, col = "black", trapnos = FALSE,
                         showcapt = length(dets) == 1 && dets != "all",
                         xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, partition)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
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
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  sigmatoa <- coefs["sigmatoa"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
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
      plot.other.contours(simpledens, toadens, plot.simple, plot.extra, D, mask)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat, col, ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}

#' @rdname contours
#' @method contours ss
#' @S3method contours ss
contours.ss <- function(fit, dets = "all", add = FALSE, heat = FALSE,
                        col = "black", trapnos = FALSE,
                        showcapt = length(dets) == 1 && dets != "all",
                        xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
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
  D <- coefs["D"]
  ssb0 <- coefs["ssb0"]
  ssb1 <- coefs["ssb1"]
  sigmass <- coefs["sigmass"]
  muss <- ssb0 + ssb1*dist
  allnonprobs <- pnorm(cutoff, muss, sigmass)
  for (i in dets){
    ssdens <- logdens.ss(allcapt, allsscapt, allnonprobs, ntraps,
                         muss, sigmass, i)
    maskdens <- exp(ssdens)*D
    maskdens <- maskdens/sum(maskdens)
    plot.main.contour(maskdens, mask, xlim, ylim, heat, col, ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}

#' @rdname contours
#' @method contours ang
#' @S3method contours ang
contours.ang <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                         heat = FALSE, col = "black", trapnos = FALSE,
                         showcapt = length(dets) == 1 && dets != "all",
                         xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, partition)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
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
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  kappa <- coefs["kappa"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    angdens <- logdens.ang(allangcapt, allcapt, ang, kappa, i)
    maskdens <- exp(simpledens + angdens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(simpledens, angdens, plot.simple, plot.extra, D, mask)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat, col, ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}

#' @rdname contours
#' @method contours disttc
#' @S3method contours disttc
contours.disttc <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                            heat = FALSE, col = "black", trapnos = FALSE,
                            showcapt = length(dets) == 1 && dets != "all",
                            xlim = NULL, ylim = NULL, ...){
  data <- fit$data
  n <- data$n
  updated.arguments <- warning.contours(n, dets, add, heat, showcapt, partition)
  dets <- updated.arguments$dets
  showcapt <- updated.arguments$showcapt
  partition <- updated.arguments$partition
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
  D <- coefs["D"]
  g01 <- coefs["g01"]
  sigma1 <- coefs["sigma1"]
  g02 <- coefs["g02"]
  sigma2 <- coefs["sigma2"]
  sigma <- coefs["sigma"]
  alpha <- coefs["alpha"]
  allprobs <- matrix(0, nrow = nrow(dist), ncol = ncol(dist))
  allprobs[1, ] <- g01*exp(-dist[1, ]^2/(2*sigma1^2))
  allprobs[2, ] <- g02*exp(-dist[2, ]^2/(2*sigma2^2))
  for (i in dets){
    simpledens <- logdens.simple(allcapt, allprobs, ntraps, i)
    distdens <- logdens.disttc(alldistcapt, allcapt, dist, alpha, i)
    maskdens <- exp(simpledens + distdens)*D
    maskdens <- maskdens/sum(maskdens)
    if (plot.part){
      plot.other.contours(simpledens, distdens, plot.simple, plot.extra, D, mask)
    }
    plot.main.contour(maskdens, mask, xlim, ylim, heat, col, ...)
  }
  plot.traps(traps, allcapt, i, heat, trapnos, showcapt)
}


## Checks inputs and returns altered argument values.
warning.contours <- function(n, dets, add, heat, showcapt, partition = NULL){
  if (length(dets) == 1 && dets == "all"){
    dets <- 1:n
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
  if (!is.null(partition)){
    if (heat & !(partition == FALSE | partition == "none")){
      warning("Setting partition to FALSE as heat is TRUE")
      partition <- FALSE
    }
    if (length(dets) > 1 & !(partition == FALSE | partition == "none")){
      warning("Setting partition to FALSE as length(dets) > 1")
      partition <- FALSE
    }
  }
  list(dets = dets, showcapt = showcapt, partition = partition)
}

## Sets up indicators for simple and extra contour plotting.
check.partition <- function(partition){
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
  if (require(TeachingDemos)){
    op <- TeachingDemos::squishplot(xlim, ylim, 1)
    plot(mask, type = "n", xlim = xlim, ylim = ylim)
    par(op)
  } else {
    plot(mask, type = "n", xlim = xlim, ylim = ylim)
    warning("Make package 'TeachingDemos' available to ensure an aspect ratio of 1.")
  }
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

logdens.disttc <- function(alldistcapt, allcapt, dist, alpha, i){
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
plot.main.contour <- function(maskdens, mask, xlim, ylim, heat, col, ...){
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
          xlim = xlim, ylim = ylim)
    box()
  } else {
    contour(x = unique.x, y = unique.y, z = z, add = TRUE, col = col, ...)
  }
}

## Plots the extra contours (i.e., when partition != FALSE).
plot.other.contours <- function(detdens, otherdens, plot.simple,
                                plot.extra, D, mask){
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
  z1col <- rgb(0, 1, 0, 0.4)
  z2col <- rgb(0, 0, 1, 0.4)
  if (plot.simple){
      contour(x = unique.x, y = unique.y, z = z1, add = TRUE, col = z1col)
  }
  if (plot.extra){
      contour(x = unique.x, y = unique.y, z = z2, add = TRUE, col = z2col)
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
