#' Plotting estimated animal location densities.
#'
#' Plots a density of individual animals' locations from a fit returned by
#' \code{admbsecr()}.
#'
#' @param fit a fitted model returned by \code{admbsecr()}.
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
#' @inheritParams graphics::plot.window
#' @method contours simple
#' @S3method contours simple
contours.simple <- function(fit, dets = "all", add = FALSE, heat = FALSE,
                            col = "black", trapnos = FALSE,
                            showcapt = length(dets) == 1 && dets != "all",
                            xlim = NULL, ylim = NULL, ...){
  if (heat & add){
     warning("Setting add to FALSE as heat is TRUE")
  }
  if (heat & (length(dets) > 1)){
     stop("Only one animal can be plotted when heat is TRUE")
  }
  if (length(dets) == 1 && dets == "all"){
    dets <- 1:fit$data$n
  }
  if (length(dets) > 1 & showcapt){
     warning("Setting showcapt to FALSE as length(dets) > 1")
     showcapt <- FALSE
  }
  data <- fit$data
  mask <- fit$mask
  allcapt <- data$capt
  traps <- fit$traps
  dist <- data$dist
  n <- data$n
  ntraps <- data$ntraps
  coefs <- coef(fit)
  x <- mask[, 1]
  y <- mask[, 2]
  if (is.null(xlim)) xlim <- range(x)
  if (is.null(ylim)) ylim <- range(y)
  if (!add & !heat){
    if (require(TeachingDemos)){
      op <- TeachingDemos::squishplot(xlim, ylim,
                                        diff(ylim)/diff(xlim))
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
      par(op)
    } else {
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
      warning("Make package 'TeachingDemos' available to ensure an aspect ratio of 1.")
    }
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  for (i in dets){
    capt <- allcapt[i, ]
    probs <- allprobs
    for (j in 1:ntraps){
      if (capt[j] == 0) probs[j, ] <- 1 - probs[j, ]
    }
    maskprobs <- exp(apply(log(probs), 2, sum))*D
    maskprobs <- maskprobs/sum(maskprobs)
    uniquex <- sort(unique(x))
    uniquey <- sort(unique(y))
    z <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
    for (j in 1:length(maskprobs)){
      xind <- which(uniquex == x[j])
      yind <- which(uniquey == y[j])
      z[xind, yind] <- maskprobs[j]
    }
    if (heat){
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",
            xlim = xlim, ylim = ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col, ...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol)
  }
  if(showcapt) {
    points(traps[which(fit$data$capt[dets, ] == 1), , drop = FALSE], cex = 2,
           lwd = 2, col = trapcol)
  }
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
  if (length(dets) == 1 && dets == "all"){
    dets <- 1:fit$data$n
  }
  if (heat & add){
     warning("Setting add to FALSE as heat is TRUE")
  }
  if (heat & !(partition == FALSE | partition == "none")){
     warning("Setting partition to FALSE as heat is TRUE")
  }
  if (heat & (length(dets) > 1)){
     stop("Only one animal can be plotted when heat is TRUE")
  }
  if (!(partition == FALSE | partition == "none") & length(dets) > 1){
    warning("Setting partition to FALSE as length(dets) > 1")
    partition <- FALSE
  }
  if (length(dets) > 1 & showcapt){
     warning("Setting showcapt to FALSE as length(dets) > 1")
     showcapt <- FALSE
  }
  if (is.logical(partition)){
    plot.simple <- partition
    plot.extra <- partition
    plot.part <- partition
  } else {
    plot.part <- TRUE
    if (partition == "all"){
      plot.simple <- TRUE
      plot.extra <- TRUE
    } else if (partition == "none"){
      plot.simple <- FALSE
      plot.extra <- FALSE
      plot.part <- FALSE
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
  data <- fit$data
  mask <- fit$mask
  allcapt <- data$capt
  alltoacapt <- data$toacapt
  traps <- fit$traps
  dist <- data$dist
  n <- data$n
  ntraps <- data$ntraps
  coefs <- coef(fit)
  x <- mask[, 1]
  y <- mask[, 2]
  if (is.null(xlim)){
    xlim <- range(x)
  }
  if (is.null(ylim)){
    ylim <- range(y)
  }
  if (!add & !heat){
    if (require(TeachingDemos)){
      op <- TeachingDemos::squishplot(xlim, ylim,
                                        diff(ylim)/diff(xlim))
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
      par(op)
    } else {
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
    }
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  sigmatoa <- coefs["sigmatoa"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  times <- dist/330
  for (i in dets){
    capt <- allcapt[i, ]
    probs <- allprobs
    for (j in 1:ntraps){
      if (capt[j] == 0){
        probs[j, ] <- 1 - probs[j, ]
      }
    }
    ## Can only incorporate TOA part if more than one detection.
    if (sum(capt) > 1){
      toacapt <- alltoacapt[i, ]
      dettraps <- which(capt == 1)
      ## Getting mask point sound travel times from traps at which a
      ## detection was made.
      dettimes <- times[dettraps, ]
      toacapt <- toacapt[dettraps]
      ## Expected sound generation times.
      esttimes <- toacapt - dettimes
      ## Calculating TOA sum of squares.
      ssqtoa <- apply(esttimes, 2, function(x) sum((x - mean(x))^2))
      toadens <- (1 - sum(capt))*log(sigmatoa^2) - ssqtoa/(2*sigmatoa^2)
    } else {
      if (partition){
        warning("Setting partition to FALSE; no TOA information.")
        partition <- FALSE
      }
      toadens <- 0
    }
    maskprobs <- exp(apply(log(probs), 2, sum) + toadens)*D
    maskprobs <- maskprobs/sum(maskprobs)
    uniquex <- sort(unique(x))
    uniquey <- sort(unique(y))
    z <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
    for (j in 1:length(maskprobs)){
      xind <- which(uniquex == x[j])
      yind <- which(uniquey == y[j])
      z[xind, yind] <- maskprobs[j]
    }
    if (plot.part){
      secrprobs <- exp(apply(log(probs), 2, sum))*D
      secrprobs <- secrprobs/sum(secrprobs)
      toaprobs <- exp(toadens)*D
      toaprobs <- toaprobs/sum(toaprobs)
      z1 <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
      z2 <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
      for (j in 1:length(maskprobs)){
        xind <- which(uniquex == x[j])
        yind <- which(uniquey == y[j])
        z1[xind, yind] <- secrprobs[j]
        z2[xind, yind] <- toaprobs[j]
      }
      z1col <- rgb(0, 1, 0, 0.4)
      z2col <- rgb(0, 0, 1, 0.4)
      if (plot.simple){
        contour(x = uniquex, y = uniquey, z = z1, add = TRUE, col = z1col)
      }
      if (plot.extra){
        contour(x = uniquex, y = uniquey, z = z2, add = TRUE, col = z2col)
      }
    }
    if (heat){
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",
            xlim = xlim, ylim = ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col, ...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol)
  }
  if(showcapt) {
    points(traps[which(fit$data$capt[dets, ] == 1), , drop = FALSE], cex = 2,
           lwd = 2, col = trapcol)
  }
}

#' @rdname contours
#' @method contours ss
#' @S3method contours ss
contours.ss <- function(fit, dets = "all", add = FALSE, heat = FALSE,
                        col = "black", trapnos = FALSE,
                        showcapt = length(dets) == 1 && dets != "all",
                        xlim = NULL, ylim = NULL, ...){
  if (length(dets) == 1 && dets == "all"){
    dets <- 1:fit$data$n
  }
  if (heat & add){
     warning("Setting add to FALSE as heat is TRUE")
  }
  if (heat & (length(dets) > 1)){
     stop("Only one animal can be plotted when heat is TRUE")
  }
  if (length(dets) > 1 & showcapt){
     warning("Setting showcapt to FALSE as length(dets) > 1")
     showcapt <- FALSE
  }
  data <- fit$data
  mask <- fit$mask
  allcapt <- data$capt
  allsscapt <- data$sscapt
  traps <- fit$traps
  dist <- data$dist
  n <- data$n
  ntraps <- data$ntraps
  nmask <- data$nmask
  cutoff <- fit$data$c
  coefs <- coef(fit)
  x <- mask[, 1]
  y <- mask[, 2]
  if (is.null(xlim)){
    xlim <- range(x)
  }
  if (is.null(ylim)){
    ylim <- range(y)
  }
  if (!add & !heat){
    if (require(TeachingDemos)){
      op <- TeachingDemos::squishplot(xlim, ylim,
                                        diff(ylim)/diff(xlim))
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
      par(op)
    } else {
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
    }
  }
  D <- coefs["D"]
  ssb0 <- coefs["ssb0"]
  ssb1 <- coefs["ssb1"]
  sigmass <- coefs["sigmass"]
  muss <- ssb0 + ssb1*dist
  allnonprobs <- pnorm(cutoff, muss, sigmass)
  for (i in dets){
    capt <- allcapt[i, ]
    sscapt <- allsscapt[i, ]
    probs <- matrix(0, nrow = ntraps, ncol = nmask)
    for (j in 1:ntraps){
      if (capt[j] == 1){
        probs[j, ] <- dnorm(sscapt[j], muss[j, ], sigmass)
      } else {
        probs[j, ] <- allnonprobs[j, ]
      }
    }
    maskprobs <- exp(apply(log(probs), 2, sum))*D
    maskprobs <- maskprobs/sum(maskprobs)
    uniquex <- sort(unique(x))
    uniquey <- sort(unique(y))
    z <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
    for (j in 1:length(maskprobs)){
      xind <- which(uniquex == x[j])
      yind <- which(uniquey == y[j])
      z[xind, yind] <- maskprobs[j]
    }
    if (heat){
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",
            xlim = xlim, ylim = ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col, ...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol)
  }
  if(showcapt) {
    points(traps[which(fit$data$capt[dets, ] == 1), , drop = FALSE], cex = 2,
           lwd = 2, col = trapcol)
  }
}

#' @rdname contours
#' @method contours ang
#' @S3method contours ang
contours.ang <- function(fit, dets = "all", add = FALSE, partition = FALSE,
                         heat = FALSE, col = "black", trapnos = FALSE,
                         showcapt = length(dets) == 1 && dets != "all",
                         xlim = NULL, ylim = NULL, ...){
  if (length(dets) == 1 && dets == "all"){
    dets <- 1:fit$data$n
  }
  if (heat & add){
     warning("Setting add to FALSE as heat is TRUE")
  }
  if (heat & !(partition == FALSE | partition == "none")){
     warning("Setting partition to FALSE as heat is TRUE")
  }
  if (heat & (length(dets) > 1)){
     stop("Only one animal can be plotted when heat is TRUE")
  }
  if (!(partition == FALSE | partition == "none") & length(dets) > 1){
    warning("Setting partition to FALSE as length(dets) > 1")
    partition <- FALSE
  }
  if (length(dets) > 1 & showcapt){
     warning("Setting showcapt to FALSE as length(dets) > 1")
     showcapt <- FALSE
  }
  if (is.logical(partition)){
    plot.simple <- partition
    plot.extra <- partition
    plot.part <- partition
  } else {
    plot.part <- TRUE
    if (partition == "all"){
      plot.simple <- TRUE
      plot.extra <- TRUE
    } else if (partition == "none"){
      plot.simple <- FALSE
      plot.extra <- FALSE
      plot.part <- FALSE
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
  data <- fit$data
  mask <- fit$mask
  allcapt <- data$capt
  allangcapt <- data$angcapt
  traps <- fit$traps
  dist <- data$dist
  ang <- data$ang
  n <- data$n
  ntraps <- data$ntraps
  coefs <- coef(fit)
  x <- mask[, 1]
  y <- mask[, 2]
  if (is.null(xlim)){
    xlim <- range(x)
  }
  if (is.null(ylim)){
    ylim <- range(y)
  }
  if (!add & !heat){
    if (require(TeachingDemos)){
      op <- TeachingDemos::squishplot(xlim, ylim,
                                        diff(ylim)/diff(xlim))
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
      par(op)
    } else {
      plot(mask, type = "n", xlim = xlim, ylim = ylim)
    }
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  kappa <- coefs["kappa"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  for (i in dets){
    capt <- allcapt[i, ]
    probs <- allprobs
    for (j in 1:ntraps){
      if (capt[j] == 0){
        probs[j, ] <- 1 - probs[j, ]
      }
    }
    angcapt <- allangcapt[i, ]
    dettraps <- which(capt == 1)
    angdens <- dvm(angcapt[dettraps], ang[dettraps, ], kappa)
    if (length(dettraps) == 1){
      dim(angdens) <- c(1, length(x))
    }
    angdens <- apply(log(angdens), 2, sum)
    maskprobs <- exp(apply(log(probs), 2, sum) + angdens)*D
    maskprobs <- maskprobs/sum(maskprobs)
    uniquex <- sort(unique(x))
    uniquey <- sort(unique(y))
    z <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
    for (j in 1:length(maskprobs)){
      xind <- which(uniquex == x[j])
      yind <- which(uniquey == y[j])
      z[xind, yind] <- maskprobs[j]
    }
    if (plot.part){
      secrprobs <- exp(apply(log(probs), 2, sum))*D
      secrprobs <- secrprobs/sum(secrprobs)
      angprobs <- exp(angdens)*D
      angprobs <- angprobs/sum(angprobs)
      z1 <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
      z2 <- matrix(NA, nrow = length(uniquex), ncol = length(uniquey))
      for (j in 1:length(maskprobs)){
        xind <- which(uniquex == x[j])
        yind <- which(uniquey == y[j])
        z1[xind, yind] <- secrprobs[j]
        z2[xind, yind] <- angprobs[j]
      }
      z1col <- rgb(0, 1, 0, 0.4)
      z2col <- rgb(0, 0, 1, 0.4)
      if (plot.simple){
        contour(x = uniquex, y = uniquey, z = z1, add = TRUE, col = z1col)
      }
      if (plot.extra){
        contour(x = uniquex, y = uniquey, z = z2, add = TRUE, col = z2col)
      }
    }
    if (heat){
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",
            xlim = xlim, ylim = ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col, ...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol)
  }
  if(showcapt) {
    points(traps[which(fit$data$capt[dets, ] == 1), , drop = FALSE], cex = 2,
           lwd = 2, col = trapcol)
  }
}
