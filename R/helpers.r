library(TeachingDemos) # needed for squishplot() to set aspect ratio of plots.
# NOTE: squishplot() sets margins that can screw up subsequent plots. Really need to unset
#       margins after using squishplot(). Not yet done.
library(CircStats) # needed for dvm()

contours <- function(fit, det = "all", add = FALSE, partition = FALSE,
                         heat = FALSE, col = "black", trapnos = FALSE,showcapt=TRUE,
                         xlim=NULL,ylim=NULL,nosimple=FALSE,...){
  #----------------------------------------------------------------------------------
  # xlim, ylim: x- and y- ranges for plotting
  # ...       : additional arguments to pass to function contour()
  # nosimple  : suppresses simple SECR method contours when partition=TRUE
  # showcapt  : if TRUE and length(det)==1, circles are drawn around detectors on which
  #             the detection was made.
  #----------------------------------------------------------------------------------
  method <- fit$method
  if (length(det) == 1 && det == "all"){
    det <- 1:fit$data$n
  }
  if (heat & add){
    warning("Setting add to FALSE as heat is TRUE")
  }
  if (heat & partition){
    warning("Setting partition to FALSE as heat is TRUE")
  }
  if (heat & length(det) > 1){
    stop("Only one animal can be plotted when heat is TRUE")
  }
  if (method == "simple"){
    if (partition){
      warning("When method is \"simple\" partition = TRUE is not relevant.")
    }
    if (nosimple){
      warning("When method is \"simple\" nosimple = TRUE is not relevant.")
    }
    contours.simple(fit, det, add, heat, col, trapnos,showcapt,xlim,ylim,...)
  } else if (method == "toa"){
    contours.toa(fit, det, add, partition, heat, col, trapnos,showcapt,xlim,ylim,nosimple,...)
  } else if (method == "ss"){
    if (partition){
      warning("When method is \"ss\" partition = TRUE is not relevant.")
    }
    if (nosimple){
      warning("When method is \"simple\" nosimple = TRUE is not relevant.")
    }
    contours.ss(fit, det, add, heat, col, trapnos,showcapt,xlim,ylim,...)
  } else if (method == "ang"){
    contours.ang(fit, det, add, partition, heat, col, trapnos,showcapt,xlim,ylim,nosimple,...)
  } else {
    stop(paste("This function does not work with admbsecr fits of method ",
               "\"", method, "\"", sep = ""))
  }
}

contours.simple <- function(fit, det, add, heat, col, trapnos,showcapt=TRUE,xlim,ylim,...){
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
  if(is.null(xlim)) xlim=range(x)
  if(is.null(ylim)) ylim=range(x)
  if (!add){
    squishplot(xlim,ylim,diff(ylim)/diff(xlim))
    plot(mask, type = "n",xlim=xlim,ylim=ylim)
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  for (i in det){
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
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",xlim=xlim,ylim=ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col,xlim=xlim,ylim=ylim,...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol,lwd=2)
  }
  if(showcapt) {
    if(length(det)!=1) {
      warning("Not showing capture locations because length(det)>1.")
    } else {
      points(traps[which(fit$data$capt[det,]==1),,drop=FALSE],cex=2,lwd=2,col=trapcol)
    }
  }
}

contours.toa <- function(fit, det, add, partition, heat, col, trapnos,showcapt=TRUE,xlim,ylim,nosimple,...){
  if (partition & length(det) > 1){
    warning("Setting partition to FALSE as length(det) > 1")
    partition <- FALSE
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
  if(is.null(xlim)) xlim=range(x)
  if(is.null(ylim)) ylim=range(x)
  if (!add){
    squishplot(xlim,ylim,diff(ylim)/diff(xlim))
    plot(mask, type = "n",xlim=xlim,ylim=ylim)
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  sigmatoa <- coefs["sigmatoa"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  times <- dist/330
  for (i in det){
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
    if (partition){
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
      z1col <- "blue" # rgb(0, 1, 0, 0.4)
      z2col <- "black" # rgb(0, 0, 1, 0.4)
      if(!nosimple) contour(x = uniquex, y = uniquey, z = z1, add = TRUE, col = z1col,xlim=xlim,ylim=ylim,lty=2)
      contour(x = uniquex, y = uniquey, z = z2, add = TRUE, col = z2col,xlim=xlim,ylim=ylim,lty=3)
    }
    if (heat){
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",xlim=xlim,ylim=ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col,xlim=xlim,ylim=ylim,...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol,lwd=2)
  }
  if(showcapt) {
    if(length(det)!=1) {
      warning("Not showing capture locations because length(det)>1.")
    } else {
      points(traps[which(fit$data$capt[det,]==1),,drop=FALSE],cex=2,lwd=2,col=trapcol)
    }
  }
}


contours.ss <- function(fit, det, add, heat, col, trapnos,showcapt=TRUE,xlim,ylim,...){
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
  if(is.null(xlim)) xlim=range(x)
  if(is.null(ylim)) ylim=range(x)
  if (!add){
    squishplot(xlim,ylim,diff(ylim)/diff(xlim))
    plot(mask, type = "n",xlim=xlim,ylim=ylim)
  }
  D <- coefs["D"]
  ssb0 <- coefs["ssb0"]
  ssb1 <- coefs["ssb1"]
  sigmass <- coefs["sigmass"]
  muss <- ssb0 + ssb1*dist
  allnonprobs <- pnorm(cutoff, muss, sigmass)
  for (i in det){
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
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",xlim=xlim,ylim=ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col,xlim=xlim,ylim=ylim,...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol,lwd=2)
  }
  if(showcapt) {
    if(length(det)!=1) {
      warning("Not showing capture locations because length(det)>1.")
    } else {
      points(traps[which(fit$data$capt[det,]==1),,drop=FALSE],cex=2,lwd=2,col=trapcol)
    }
  }  
}

contours.ang <- function(fit, det, add, partition, heat, col, trapnos,showcapt=TRUE,xlim,ylim,nosimple,...){
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
  if(is.null(xlim)) xlim=range(x)
  if(is.null(ylim)) ylim=range(x)
  if (!add){
    squishplot(xlim,ylim,diff(ylim)/diff(xlim))
    plot(mask, type = "n",xlim=xlim,ylim=ylim)
  }
  D <- coefs["D"]
  g0 <- coefs["g0"]
  sigma <- coefs["sigma"]
  kappa <- coefs["kappa"]
  allprobs <- g0*exp(-dist^2/(2*sigma^2))
  for (i in det){
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
    if (partition){
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
      z1col <- "blue" # rgb(0, 1, 0, 0.4)
      z2col <- "black" # rgb(0, 0, 1, 0.4)
      if(!nosimple) contour(x = uniquex, y = uniquey, z = z1, add = TRUE, col = z1col,xlim=xlim,ylim=ylim,lty=2)
      contour(x = uniquex, y = uniquey, z = z2, add = TRUE, col = z2col,xlim=xlim,ylim=ylim,lty=3)
    }
    if (heat){
      image(x = uniquex, y = uniquey, z = z, xlab = "x", ylab = "y",xlim=xlim,ylim=ylim)
      box()
      trapcol <- "black"
    } else {
      contour(x = uniquex, y = uniquey, z = z, add = TRUE, col = col,xlim=xlim,ylim=ylim,...)
      trapcol <- "red"
    }
  }
  if (trapnos){
    text(traps, labels = 1:ntraps, col = trapcol)
  } else {
    points(traps, pch = 4, col = trapcol,lwd=2)
  }
  if(showcapt) {
    if(length(det)!=1) {
      warning("Not showing capture locations because length(det)>1.")
    } else {
      points(traps[which(fit$data$capt[det,]==1),,drop=FALSE],cex=2,lwd=2,col=trapcol)
    }
  }  
}

