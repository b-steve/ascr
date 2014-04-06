#' Simulating SECR data
#'
#' Simulates SECR capture histories and associated additional
#' information in the correct format for use with the function
#' \link{admbsecr}. If \code{fit} is provided then no other arguments
#' are required. Otherwise, at least \code{traps}, \code{mask}, and
#' \code{pars} are needed.
#'
#' See documentation for the function \link{admbsecr} for information
#' on the parameters corresponding to the different detection
#' functions, and to different types of additional information.
#'
#' Simulated call frequencies are not always integers, e.g. when
#' \code{freq.dist} is \code{"norm"}, or when \code{freq.dist} is
#' \code{"edf"} and the call frequencies used to fit the model
#' \code{fit} are not all integers. In this case, if \code{freq.dist}
#' is \code{"edf"}, then simulated call frequencies are rounded at
#' random as follows: Let \eqn{x} be the fraction part of the number,
#' then the call frequency is rounded up with probability \eqn{x} and
#' rounded down with probability \eqn{1 - x}. For example, a value of
#' 8.1 will be rounded to 9 with probability 0.1, and rounded to 8
#' with probability 0.9.
#'
#' @param fit A fitted \code{admbsecr} model object which provides the
#' additional information types, detection function, and parameter
#' values from which to generate capture histories.
#' @param infotypes A character vector indicating the type(s) of
#' additional information to be simulated. Elements can be a subset of
#' \code{"bearing"}, \code{"dist"}, \code{"ss"}, \code{"toa"}, and
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
#' @param freq.dist A character string, either \code{"edf"} or
#' \code{"norm"}, which specifies how the distribution function of the
#' call frequencies should be estimated. If \code{"edf"}, then the
#' distribution of call frequencies is estimated using the empirical
#' distribution function. If \code{"norm"}, then a normal distribution
#' is fitted to the call frequencies using the sample mean and
#' variance. See 'Details' below for information on how call
#' frequencies are rounded.
#' @param test.detfn Logical value, if \code{TRUE}, tests detection
#' function to aid debugging.
#' @inheritParams admbsecr
#'
#' @return A list with named components, each corresponding to a data
#' type listed in \code{infotypes}. Each component is a matrix where
#' each row corresponds to each detected individual, and each column
#' corresponds to a trap (or detector). The elements in the matrix
#' indicate detection, and provide simulated values of the additional
#' information requested. This object can be used as the \code{capt}
#' argument for the function \link{admbsecr}.
#' 
#' @examples
#' ## Simulating based on model fits.
#' simple.capt <- sim.capt(simple.hn.fit)
#' bearing.capt <- sim.capt(bearing.hn.fit)
#' ## Simulating from provided parameter values.
#' new.capt <- sim.capt(traps = example.traps, mask = example.mask, infotypes = c("bearing", "dist"), detfn = "hr",
#'                      pars = list(D = 2500, g0 = 0.9, sigma = 3, z = 2, kappa = 50, alpha = 10))
#' 
#' @export
sim.capt <- function(fit = NULL, traps = NULL, mask = NULL,
                     infotypes = character(0), detfn = "hn",
                     pars = NULL, ss.link = "identity",
                     cutoff = NULL, call.freqs = NULL, freq.dist = "edf",
                     sound.speed = 330, test.detfn = FALSE){
    ## Some error checking.
    if (any(infotypes == "ss")){
        stop("Signal strength information is simulated by setting argument 'detfn' to \"ss\".")
    }
    if (!missing(cutoff) & detfn != "ss"){
        warning("The argument 'cutoff' is being ignored, as 'detfn' is not \"ss\".")
    }
    if (!missing(ss.link) & detfn != "ss"){
        warning("The argument 'ss.link' is being ignored, as 'detfn' is not \"ss\".")
    }
    ## Grabbing values from fit if required.
    if (!is.null(fit)){
        traps <- get.traps(fit)
        mask <- get.mask(fit)
        infotypes <- fit$infotypes
        detfn <- fit$args$detfn
        pars <- get.par(fit, "fitted", as.list = TRUE)
        ss.link <- fit$args$ss.link
        cutoff <- fit$args$cutoff
        call.freqs <- fit$args$call.freqs
        sound.speed <- fit$args$sound.speed
    }
    ## Setting up logical indicators for additional information types.
    supp.types <- c("bearing", "dist", "ss", "toa", "mrds")
    sim.types <- supp.types %in% infotypes
    names(sim.types) <- supp.types
    sim.bearings <- sim.types["bearing"]
    sim.dists <- sim.types["dist"]
    sim.toas <- sim.types["toa"]
    sim.mrds <- sim.types["mrds"]
    sim.ss <- ifelse(detfn == "ss", TRUE, FALSE)
    ## Working out required parameters.
    suppar.names <- c("kappa", "alpha", "sigma.toa")[sim.types[c("bearing", "dist", "toa")]]
    if (sim.ss){
        if (ss.link == "identity"){
            detfn <- "ss"
        } else if (ss.link == "log"){
            detfn <- "log.ss"
        } else {
            stop("The argument 'ss.link' must be either \"identity\" or \"log\"")
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
    ## Grabbing detection function parameters.
    detpars <- pars[detpar.names]
    ## Specifies the area in which animal locations can be generated.
    core <- data.frame(x = range(mask[, 1]), y = range(mask[, 2]))
    ## Simulating population.
    if (is.null(call.freqs)){
        popn <- as.matrix(sim.popn(D = pars$D, core = core, buffer = 0))
    } else {
        D <- pars$D/mean(call.freqs)
        popn <- as.matrix(sim.popn(D = D, core = core, buffer = 0))
        n.a <- nrow(popn)
        if (freq.dist == "edf"){
            if (length(call.freqs) == 1){
                freqs <- rep(call.freqs, n.a)
            } else {
                freqs <- sample(call.freqs, size = n.a, replace = TRUE)
            }
        } else if (freq.dist == "norm"){
            if (diff(range(call.freqs)) == 0){
                freqs <- rep(unique(call.freqs), n.a)
            } else {
                freqs <- rnorm(n.a, mean(call.freqs), sd(call.freqs))
            }
        } else {
            stop("The argument 'freq.dist' must be either \"edf\" or \"norm\"")
        }
        ## Rounding frequencies up and down at random, depending
        ## on which integer is closer.
        which.integers <- floor(freqs) == freqs
        for (i in (1:n.a)[!which.integers]){
            prob <- freqs[i] - floor(freqs[i])
            freqs[i] <- floor(freqs[i]) + rbinom(1, 1, prob)
        }
        popn <- popn[rep(1:n.a, times = freqs), ]
    }
    n.popn <- nrow(popn)
    if (n.popn == 0) stop("No animals in population.")
    ## Calculating distances.
    dists <- distances(popn, traps)
    n.traps <- nrow(traps)
    ## Calculating detection probabilities and simulating captures.
    if (!sim.ss){
        det.probs <- calc.detfn(dists, detfn, detpars, ss.link)
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
        lines(xx, calc.detfn(xx, detfn, detpars), col = "blue")
    }
    ## Total number of detections.
    n.dets <- sum(bin.capt)
    ## Locations of captured individuals.
    capt.popn <- popn[captures, ]
    ## Capture distances.
    capt.dists <- dists[captures, ]
    ## Simulating additional information.
    if (sim.bearings){
        bearings <- t(bearings(traps, capt.popn))
        bearing.capt <- matrix(0, nrow = nrow(bin.capt),
                           ncol = ncol(bin.capt))
        bearing.capt[bin.capt == 1] <- (bearings[bin.capt == 1] +
                     rvm(n.dets, mean = 0, k = pars$kappa)) %% (2*pi)
        out$bearing <- bearing.capt
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

