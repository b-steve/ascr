## Returns capture trap numbers.
trapvec <- function(capthist){
  x <- apply(capthist, 3, function(x) sum(x > 0))
  rep(1:length(x), times = x)
}

## Returns capture animal ID numbers.
animalIDvec <- function(capthist){
  x <- c(apply(capthist, 3, function(x) which(x > 0)), recursive = TRUE)
  names(x) <- NULL
  as.character(x)
}

#' Assigning ID numbers to sounds.
#'
#' Identifies recaptures and assigns ID numbers to sounds recorded for
#' an SECR model.
#'
#' Detected sounds are assumed to come from the same animal if times
#' of arrival at different microphones are closer together than the
#' time it would take for sound to travel between these microphones.
#'
#' @param mics a matrix containing the coordinates of trap locations.
#' @param clicks a data frame containing (at least): (i) \code{$tim$},
#' the precise time of arrival of the received sound, and (ii)
#' \code{$trap} the trap at which the sound was recorded.
#' @param dt a \code{K} by \code{K} matrix (where \code{K} is the
#' number of traps) containing the time taken for sound to travel
#' between each pair of traps.
#' @return A data frame. Specifically, the \code{clicks} dataframe,
#' now with a new variable, \code{ID}.
#' @author David Borchers, Ben Stevenson
#' @export
make.acoustic.captures <- function(mics, clicks, dt){
  K <- dim(mics)[1]
  captures <- clicks
  ct <- rep(-Inf, K)
  ID <- 1
  ct[clicks$trap[1]] <- clicks$tim[1]
  new <- FALSE
  nclicks <- length(clicks$tim)
  for (i in 2:nclicks){
    if (ct[clicks$trap[i]] > -Inf){
      nd <- length(which(ct > -Inf))
      captures$ID[(i - nd):(i - 1)] <- ID
      ct <- rep(-Inf, K)
      ct[clicks$trap[i]] <- clicks$tim[i]
      ID <- ID + 1
      if(i == nclicks) captures$ID[i] <- ID
    }
    else {
      ct[clicks$trap[i]] <- clicks$tim[i]
      ctset <- which(ct > -Inf)
      dts <- dt[ctset, clicks$trap[i]]
      cts <- -(ct[ctset] - clicks$tim[i])
      if (any((cts - dts) > 0)) new <- TRUE
      if (new) {
        nd <- length(which(ct > -Inf)) - 1
        captures$ID[(i - nd):(i - 1)] <- ID
        ct <- rep(-Inf, K)
        ct[clicks$trap[i]] <- clicks$tim[i]
        ID <- ID + 1
        new <- FALSE
        if (i == nclicks) captures$ID[i] <- ID
      } else if(i == nclicks){
        nd <- length(which(ct > -Inf))
        captures$ID[(i - nd + 1):i] <- ID
      }
    }
  }
  captures
}

## Adapted from R2admb.
read.admbsecr <- function(fn, verbose = FALSE, checkterm = TRUE){
  if (verbose)
    cat("reading output ...\n")
  parfn <- paste(fn, "par", sep = ".")
  if (!file.exists(parfn))
    stop("couldn't find parameter file ", parfn)
  L <- c(list(fn = fn), read_pars(fn))
  if (checkterm) {
    v <- with(L, vcov[seq(npar), seq(npar)])
    ev <- try(eigen(solve(v))$value, silent = TRUE)
    L$eratio <- if (inherits(ev, "try-error"))
      NA
    else min(ev)/max(ev)
  }
  class(L) <- "admb"
  L
}

#' Extract parameter standard errors.
#'
#' Extracts standard errors from an admbsecr fit.
#'
#' @param fit a fitted model from \code{admbsecr()}.
#' @export
stdEr <- function(fit){
    fit$se
}

## Return fixed or estimated parameter values from a model fit.
getpar <- function(fit, pars, as.list = FALSE){
    allpar.names <- c("D", fit$detpars, fit$suppars)
    if (length(pars) == 1){
        if (pars == "all"){
            pars <- allpar.names
        }
    }
    ## Error checking.
    legal.names <- pars %in% allpar.names
    if (!all(legal.names)){
        illegal.pars <- pars[!legal.names]
        if (sum(!legal.names) == 1){
            msg <- paste(illegal.pars, "is not a parameter in the model provided.")
        } else if (sum(!legal.names) == 2){
            msg <- paste(paste(illegal.pars, collapse = " and "),
                         "are not parameters in the model provided.")
        } else if (sum(!legal.names) > 2){
            n.illegal <- length(illegal.pars)
            msg <- paste(paste(illegal.pars[-n.illegal], collapse = ", "),
                         ", and", illegal.pars[n.illegal],
                         "are not parameters in the model provided.")
        }
        stop(msg)
    }
    out <- numeric(length(pars))
    names(out) <- pars
    det.index <- which(fit$detpars %in% pars)
    supp.index <- which(fit$suppars %in% pars)
    ## Logical vector indicating parameters that weren't estimated.
    fixed.pars <- fit$phases[pars] == -1
    ## Putting in fixed parameter values.
    if (sum(fixed.pars) > 0){
        out[fixed.pars] <- fit$sv[[pars[fixed.pars]]]
    }
    ## Working out parameter groups for parameters in 'pars'.
    det.index <- which(pars %in% fit$detpars)
    supp.index <- which(pars %in% fit$suppars)
    admb.pars <- pars
    admb.pars[det.index] <- paste("detpars[", which(fit$detpars %in% pars), "]",
                                  sep = "")
    admb.pars[supp.index] <- paste("suppars[", which(fit$suppars %in% pars), "]",
                                  sep = "")
    ## Putting in estimated parameter values.
    out[!fixed.pars] <- coef(fit)[admb.pars[!fixed.pars]]
    if (as.list){
        out.vec <- out
        out <- vector("list", length = length(out.vec))
        names(out) <- names(out.vec)
        names(out.vec) <- NULL
        for (i in 1:length(out.vec)){
            out[[i]] <- out.vec[i]
        }
    }
    out
}

#' Simulating SECR data
#'
#' Simulates SECR capture histories and associated additional
#' information in the correct format for use with
#' \code{\link[admbsecr]{admbsecr}}.
#'
#' If \code{fit} is provided then no other arguments are
#' required. Otherwise, at least \code{traps}, \code{mask}, and
#' \code{pars} are needed.
#'
#' @param fit A fitted \code{admbsecr} model object which provides the
#' additional information types, detection function, and parameter
#' values from which to generate capture histories.
#' @param traps A matrix with two columns. The rows provide Cartesian
#' coordiates for trap locations.
#' @param mask A matrix with two columns. The rows provide Cartesian
#' coordinates for the mask point locations.
#' @param infotypes A character vector indicating the type(s) of
#' additional information to be simulated. Elements can be a subset of
#' \code{"ang"}, \code{"dist"}, \code{"ss"}, \code{"toa"}, and
#' \code{"mrds"} (NOTE: \code{"mrds"} not yet implemented).
#' @param detfn A character string specifying the detection function
#' to be used. Options are "hn" (halfnormal), "hr" (hazard rate), "th"
#' (threshold), "lth" (log-link threshold), or "ss" (signal strength).
#' @param pars A named list. Component names are parameter names, and
#' each component is the value of the associated parameter. A value
#' for the parameter \code{D}, animal density (or call density, if it
#' an acoustic survey) must always be provided, along with values for
#' parameters associated with the chosen detection function and
#' additional information type(s).
#' @param ss.link A character string, either \code{"indentity"} or
#' \code{"log"}, which specifies the link function for the signal
#' strength detection function. Only required when \code{detfn} is
#' \code{"ss"}.
#' @param cutoff The signal strength threshold, above which sounds are
#' identified as detections. Only required when \code{detfn} is
#' \code{"ss"}.
#' @param sound.speed The speed of sound in metres per second,
#' defaults to 330 (the speed of sound in air). Only used when
#' \code{info} includes \code{"toa"}.
#' @param test.detfn Logical value, if \code{TRUE}, tests detection
#' function to aid debugging.
#' @export
sim.capt <- function(fit = NULL, traps = NULL, mask = NULL,
                     infotypes = character(0), detfn = "hn",
                     pars = NULL, ss.link = "identity",
                     cutoff = NULL, sound.speed = 330,
                     test.detfn = FALSE){
    ## Grabbing values from fit if required.
    if (!is.null(fit)){
        traps <- fit$traps
        mask <- fit$mask
        infotypes <- fit$infotypes
        detfn <- fit$detfn
        pars <- getpar(fit, "all", as.list = TRUE)
        ss.link <- fit$ss.link
        cutoff <- fit$cutoff
        sound.speed <- fit$sound.speed
    }
    ## Grabbing detection function.
    calc.detfn <- get.detfn(detfn)
    ## Setting up logical indicators for additional information types.
    supp.types <- c("ang", "dist", "ss", "toa", "mrds")
    sim.types <- supp.types %in% infotypes
    names(sim.types) <- supp.types
    sim.angs <- sim.types["ang"]
    sim.dists <- sim.types["dist"]
    sim.toas <- sim.types["toa"]
    sim.mrds <- sim.types["mrds"]
    sim.ss <- ifelse(detfn == "ss", TRUE, FALSE)
    ## Working out required parameters.
    suppar.names <- c("kappa", "alpha", "sigma.toa")[sim.types[c("ang", "dist", "toa")]]
    if (sim.ss){
        if (ss.link == "identity"){
            detfn <- "ss"
        } else if (ss.link == "log"){
            detfn <- "log.ss"
        } else {
            stop("ss.link must be either \"identity\" or \"log\"")
        }
    }
    detpar.names <- switch(detfn,
                           hn = c("g0", "sigma"),
                           hr = c("g0", "sigma", "z"),
                           th = c("shape", "scale"),
                           lth = c("shape.1", "shape.2", "scale"),
                           ss = c("b0.ss", "b1.ss", "sigma.ss"),
                           log.ss = c("b0.ss", "b1.ss", "sigma.ss"))
    par.names <- c("D", detpar.names, suppar.names)
    if (!identical(sort(par.names), sort(names(pars)))){
        msg <- paste("The following must be named components of the list 'pars': ",
                     paste(par.names, collapse = ", "), ".", sep = "")
        stop(msg)
    }
    ## Specifies the area in which animal locations can be generated.
    core <- data.frame(x = range(mask[, 1]), y = range(mask[, 2]))
    ## Simulating population.
    popn <- as.matrix(sim.popn(D = pars$D, core = core, buffer = 0))
    n.popn <- nrow(popn)
    if (n.popn == 0) stop("No animals in population.")
    ## Calculating distances.
    dists <- distances(popn, traps)
    n.traps <- nrow(traps)
    ## Calculating detection probabilities and simulating captures.
    if (!sim.ss){
        det.probs <- calc.detfn(dists, pars)
        full.bin.capt <- matrix(as.numeric(runif(n.popn*n.traps) < det.probs),
                           nrow = n.popn, ncol = n.traps)
        captures <- which(apply(full.bin.capt, 1, sum) > 0)
        bin.capt <- full.bin.capt[captures, ]
        out <- list(bincapt = bin.capt)
    } else {
        if (ss.link == "identity"){
            inv.ss.link <- identity
        } else if (ss.link == "log"){
            inv.ss.link <- exp
        } else {
            stop("Argument 'ss.link' must be \"identity\" or \"log\".")
        }
        pars$cutoff <- cutoff
        ss.mean <- inv.ss.link(pars$b0.ss - pars$b1.ss*dists)
        ss.error <- matrix(rnorm(n.popn*n.traps, mean = 0,
                                 sd = pars$sigma.ss),
                           nrow = n.popn, ncol = n.traps)
        full.ss.capt <- ss.mean + ss.error
        captures <- which(apply(full.ss.capt, 1,
                                function(x, cutoff) any(x > cutoff),
                                cutoff = cutoff))
        full.bin.capt <- ifelse(full.ss.capt > cutoff, 1, 0)
        ss.capt <- full.ss.capt[captures, ]
        bin.capt <- ifelse(ss.capt > cutoff, 1, 0)
        ss.capt[ss.capt < cutoff] <- 0
        out <- list(bincapt = bin.capt, ss = ss.capt)
    }
    ## Plot to test correct detection simulation.
    if (test.detfn){
        capt.dists <- dists[full.bin.capt == 1]
        evade.dists <- dists[full.bin.capt == 0]
        all.dists <- c(capt.dists, evade.dists)
        capt.dummy <- c(rep(1, length(capt.dists)),
                        rep(0, length(evade.dists)))
        breaks <- seq(0, max(all.dists), length.out = 100)
        mids <- breaks[-length(breaks)] + 0.5*diff(breaks)
        breaks[1] <- 0
        split.dummy <- split(capt.dummy,
                             f = cut(all.dists, breaks = breaks))
        props <- sapply(split.dummy, mean)
        plot(mids, props, type = "l", xlim = c(0, max(all.dists)),
             ylim = c(0, 1))
        xx <- seq(0, max(all.dists), length.out = 100)
        lines(xx, calc.detfn(xx, pars), col = "blue")
    }
    ## Total number of detections.
    n.dets <- sum(bin.capt)
    ## Locations of captured individuals.
    capt.popn <- popn[captures, ]
    ## Capture distances.
    capt.dists <- dists[captures, ]
    ## Simulating additional information.
    if (sim.angs){
        bearings <- t(bearings(traps, capt.popn))
        ang.capt <- matrix(0, nrow = nrow(bin.capt),
                           ncol = ncol(bin.capt))
        ang.capt[bin.capt == 1] <- (bearings[bin.capt == 1] +
                     rvm(n.dets, mean = 0, k = pars$kappa)) %% (2*pi)
        out$ang <- ang.capt
    }
    if (sim.dists){
        dist.capt <- matrix(0, nrow = nrow(bin.capt),
                            ncol = ncol(bin.capt))
        betas <- pars$alpha/capt.dists[bin.capt == 1]
        dist.capt[bin.capt == 1] <- rgamma(n.dets, shape = pars$alpha,
                      rate = betas)
        out$dist <- dist.capt
    }
    if (sim.toas){
        ## Time taken for sound to travel from source to detector.
        toa.capt <- capt.dists/sound.speed*bin.capt
        ## Adding in TOA error.
        toa.capt[bin.capt == 1] <-
            toa.capt[bin.capt == 1] + rnorm(n.dets, sd = pars$sigma.toa)
        out$toa <- toa.capt
    }
    if (sim.mrds){
        out$mrds <- capt.dists
    }
    out
}

get.detfn <- function(detfn){
    switch(detfn, hn = calc.hn, hr = calc.hr, th = calc.th,
           lth = calc.lth, ss = calc.ss, log.ss = calc.log.ss)
}

calc.hn <- function(d, pars){
    g0 <- pars$g0
    sigma <- pars$sigma
    g0*exp(-(d^2/(2*sigma^2)))
}

calc.hr <- function(d, pars){
    g0 <- pars$g0
    sigma <- pars$sigma
    z <- pars$z
    g0*(1 - exp(-((d/sigma)^-z)))
}

calc.th <- function(d, pars){
    scale <- pars$scale
    shape <- pars$shape
    0.5 - 0.5*erf(d/scale - shape)
}

calc.lth <- function(d, pars){
    scale <- pars$scale
    shape.1 <- pars$shape.1
    shape.2 <- pars$shape.2
    0.5 - 0.5*erf(shape.1 - exp(shape.2 + scale*d))
}

calc.ss <- function(d, pars){
    b0.ss <- pars$b0.ss
    b1.ss <- pars$b1.ss
    sigma.ss <- pars$sigma.ss
    cutoff <- pars$cutoff
    1 - pnorm(cutoff, mean = b0.ss - b1.ss*d, sd = sigma.ss)
}

calc.log.ss <- function(d, pars){
    b0.ss <- pars$b0.ss
    b1.ss <- pars$b1.ss
    sigma.ss <- pars$b1.ss
    cutoff <- pars$cutoff
    1 - pnorm(cutoff, mean = exp(b0.ss - b1.ss*d), sd = sigma.ss)
}

## Error function
erf <- function(x){
    2*pnorm(x*sqrt(2)) - 1
}
