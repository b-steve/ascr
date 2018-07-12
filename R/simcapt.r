#' Simulating SECR data
#'
#' Simulates SECR capture histories and associated additional
#' information in the correct format for use with the function
#' \link{fit.ascr}. If \code{fit} is provided then no other arguments
#' are required. Otherwise, at least \code{traps}, \code{mask}, and
#' \code{pars} are needed.
#'
#' See documentation for the function \link{fit.ascr} for information
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
#' If \code{cue.rates} is \code{Inf} then all simulated individuals
#' will be detected. To generate sensible capture histories then both
#' a lower and upper cutoff must be supplied in \code{ss.opts}. In
#' this case, individuals continue to emit calls until one is detected
#' above the lower cutoff by at least one microphone. This individual
#' is then incorporated into the resulting capture history only if the
#' loudest received signal strength is also above the upper cutoff.
#'
#' @param fit A fitted \code{ascr} model object which provides the
#'     additional information types, detection function, and parameter
#'     values from which to generate capture histories.
#' @param mask A matrix with two columns, providing x- and
#'     y-coordinates respectively. The extreme x- and y-coordinates
#'     define the rectangle in which individuals' locations are
#'     simulated.
#' @param infotypes A character vector indicating the type(s) of
#'     additional information to be simulated. Elements can be a
#'     subset of \code{"bearing"}, \code{"dist"}, \code{"toa"}, and
#'     \code{"mrds"} (NOTE: \code{"mrds"} not yet implemented). If
#'     signal strength information is required, specify \code{detfn}
#'     as \code{"ss"} rather than including it here.
#' @param detfn A character string specifying the detection function
#'     to be used. Options are "hn" (halfnormal), "hr" (hazard rate),
#'     "th" (threshold), "lth" (log-link threshold), or "ss" (signal
#'     strength).
#' @param popn A matrix with two columns, providing x- and
#'     y-coordinates of the animal locations.
#' @param pars A named list. Component names are parameter names, and
#'     each component is the value of the associated parameter. A
#'     value for the parameter \code{D}, animal density (or call
#'     density, if it an acoustic survey) must always be provided,
#'     along with values for parameters associated with the chosen
#'     detection function and additional information type(s).
#' @param cue.rates A vector of call frequencies from which a
#'     distribution for the number of emitted calls for each
#'     individual is fitted. If scalar, all individuals make the same
#'     number of calls. If \code{Inf}, \code{first.only} must be
#'     \code{TRUE}, and all individuals keep making calls until the
#'     first is detected.
#' @param survey.length The length of a cue-based survey.
#' @param freq.dist A character string, either \code{"edf"} or
#'     \code{"norm"}, which specifies how the distribution function of
#'     the call frequencies should be estimated. If \code{"edf"}, then
#'     the distribution of call frequencies is estimated using the
#'     empirical distribution function. If \code{"norm"}, then a
#'     normal distribution is fitted to the call frequencies using the
#'     sample mean and variance. If \code{cue.rates} is scalar then
#'     this is ignored, and all individuals make the same number of
#'     calls. See 'Details' below for information on how call
#'     frequencies are rounded.
#' @param ihd.surf A list with a vector for each session containing a
#'     value for animal density at each mask point. Only required to
#'     simulate inhomogeneous density surfaces.
#' @param test.detfn Logical value, if \code{TRUE}, tests detection
#'     function to aid debugging.
#' @param first.only Only keep the first detection for each
#'     individual.
#' @param keep.locs Logical, if \code{TRUE}, the locations of
#'     individuals in the simulated population are retained. In this
#'     case, the capture histories and auxiliary information are kept
#'     in a component \code{capt} of the returned list, locations of
#'     detected individuals are kept in a component \code{capt.locs},
#'     and locations of all individuals in the population are kept in
#'     a component \code{popn.locs}. In this case, the capture
#'     histories and auxiliary information are kept in a component
#'     \code{capt} of the returned list, and ID numbers of detected
#'     individuals are kept in a component \code{capt.ids}
#' @param keep.ids Logical, if \code{TRUE}, the ID number of detected
#'     individuals are retained.
#' @param ... Other arguments (mostly for back-compatibility).
#' @inheritParams fit.ascr
#'
#' @return A list with named components, each corresponding to a data
#'     type listed in \code{infotypes}. Each component is a matrix
#'     where each row corresponds to each detected individual, and
#'     each column corresponds to a trap (or detector). The elements
#'     in the matrix indicate detection, and provide simulated values
#'     of the additional information requested. This object can be
#'     used as the \code{capt} argument for the function
#'     \link{fit.ascr}.
#'
#' @examples
#' ## Simulating based on model fits.
#' simple.capt <- sim.capt(example$fits$simple.hn)
#' bearing.capt <- sim.capt(example$fits$bearing.hn)
#' ## Simulating from provided parameter values.
#' new.capt <- sim.capt(traps = example$traps, mask = example$mask, infotypes = c("bearing", "dist"), detfn = "hr",
#'                      pars = list(D = 2500, g0 = 0.9, sigma = 3, z = 2, kappa = 50, alpha = 10))
#'
#' @export
sim.capt <- function(fit = NULL, traps = NULL, mask = NULL, popn = NULL,
                     infotypes = character(0), detfn = "hn",
                     pars = NULL, ss.opts = NULL, cue.rates = NULL, survey.length = NULL,
                     freq.dist = "edf", sound.speed = 330, ihd.surf = NULL, test.detfn = FALSE,
                     first.only = FALSE, keep.locs = FALSE, keep.ids = FALSE, ...){
    arg.names <- names(as.list(environment()))
    extra.args <- list(...)
    if (any(names(extra.args) == "call.freqs")){
        if (!missing(cue.rates)){
            stop("The argument `cue.rates' has replaced `call.freqs'; use only the former.")
        }
        warning("The argument `call.freqs' is deprecated; please rename to `cue.rates' instead.")
        cue.rates <- extra.args[["call.freqs"]]
    }
    ## Some error checking.
    if (any(infotypes == "ss")){
        stop("Signal strength information is simulated by setting argument 'detfn' to \"ss\".")
    }
    if (!missing(ss.opts) & detfn != "ss"){
        warning("The argument 'ss.opts' is being ignored, as 'detfn' is not \"ss\".")
    }
    if (keep.ids & is.null(cue.rates)){
        warning("IDs are necessarily different as only one call is simulated from each individual.")
    }
    ## Warnings for ignored parameters. Is there a neater way of doing
    ## this? The 'missing()' function is annoying.
    if (!missing(mask) & !missing(fit)){
        warning("Argument 'mask' is being ignored as 'fit' was provided.")
    }
    if (!missing(infotypes) & !missing(fit)){
        warning("Argument 'infotypes' is being ignored as 'fit' was provided.")
    }
    if (!missing(detfn) & !missing(fit)){
        warning("Argument 'detfn' is being ignored as 'fit' was provided.")
    }
    if (!missing(traps) & !missing(fit)){
        warning("Argument 'traps' is being ignored as 'fit' was provided.")
    }
    if (!missing(pars) & !missing(fit)){
        warning("Argument 'pars' is being ignored as 'fit' was provided.")
    }
    if (!missing(ss.opts) & !missing(fit)){
        warning("Argument 'ss.opts' is being ignored as 'fit' was provided.")
    }
    if (!missing(cue.rates) & !missing(fit)){
        warning("Argument 'cue.rates' is being ignored as 'fit' was provided.")
    }
    if (!missing(sound.speed) & !missing(fit)){
        warning("Argument 'sound.speed' is being ignored as 'fit' was provided.")
    }
    if (!missing(ihd.surf) & !missing(fit)){
        warning("Argument 'ihd.surf' is being ignored as 'fit' was provided.")
    }
    ## Grabbing values from fit if required.
    if (!is.null(fit)){
        mask <- get.mask(fit)
        traps <- get.traps(fit)
        infotypes <- fit$infotypes
        detfn <- fit$args$detfn
        pars <- get.par(fit, "fitted", as.list = TRUE)
        ss.opts <- fit$args$ss.opts
        cue.rates <- fit$args$cue.rates
        survey.length <- fit$args$survey.length
        sound.speed <- fit$args$sound.speed
        ## Setting up correct arguments for a simulating from a first-call model.
        if (!is.null(ss.opts$lower.cutoff)){
            cue.rates <- Inf
            first.only <- TRUE
        }
        if (fit$fit.ihd){
            ihd.surf <- fit$D.mask
        } else {
            ihd.surf <- NULL
        }
        popn.provided <- FALSE
    } else {
        if (is.null(popn)){
            popn.provided <- FALSE
        } else {
            popn.provided <- TRUE
        }
    }
    ## Making sure mask and trap objects are lists for multisession stuff.
    full.mask <- mask
    full.traps <- traps
    if (is.list(full.mask) != is.list(full.traps)){
        stop("The 'mask' and 'traps' objects must both be lists for multisession data, or both be matrices for single-session data.")
    }
    if (is.list(full.mask)){
        if (length(full.mask) != length(full.traps)){
            stop("The 'mask' and 'traps' objects must have the same number of components.")
        }
        n.sessions <- length(full.mask)
    } else {
        n.sessions <- 1
        full.mask <- list(full.mask)
        full.traps <- list(full.traps)
    }
    ## Setting up inhomogeneous density stuff.
    if (is.null(ihd.surf)){
        sim.ihd <- FALSE
    } else {
        sim.ihd <- TRUE
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
    cutoff <- ss.opts$cutoff
    lower.cutoff <- ss.opts$lower.cutoff
    het.source <- ss.opts$het.source
    directional <- ss.opts$directional
    ss.link <- ss.opts$ss.link
    ## Sorting out inf calls stuff.
    if (any(cue.rates == Inf)){
        if (!first.only){
            warning("Setting 'first.only' to 'TRUE' as 'cue.rates' is Inf.")
            first.only <- TRUE
        }
        if (is.null(lower.cutoff)){
            #stop("A lower cutoff must be specified if individuals call until they are detected.")
        }
        inf.calls <- TRUE
    } else {
        inf.calls <- FALSE
    }
    ## Sorting out directional calling stuff.
    if (sim.ss){
        if (is.null(cutoff)){
            stop("For signal strength models, the 'cutoff' component of 'ss.opts' must be specified.")
        }
        if (is.null(directional)){
            if ("b2.ss" %in% names(pars)){
                directional <- TRUE
            } else {
                directional <- FALSE
                pars$b2.ss <- 0
            }
        } else if (directional & !("b2.ss" %in% names(pars))){
            stop("Parameter 'b2.ss' must be specified for a directional calling model.")
        } else if (!directional & "b2.ss" %in% names(pars)){
            if (pars$b2.ss != 0){
                warning("Parameter 'b2.ss' in 'pars' is being ignored as the 'directional' component of 'ss.opts' is 'FALSE'.")
                pars$b2.ss <- 0
            }
        }
        if (is.null(het.source)){
            if ("sigma.b0.ss" %in% names(pars)){
                het.source <- TRUE
            } else {
                het.source <- FALSE
                pars$sigma.b0.ss <- 0
            }
        } else if (het.source & !("sigma.b0.ss" %in% names(pars))){
            stop("Parameter 'sigma.b0.ss' must be specified for a model with heterogeneity in source strengths'.")
        } else if (!het.source & "sigma.b0.ss" %in% names(pars)){
            if (pars$sigma.b0.ss != 0){
                warning("Parameter 'sigma.b0.ss' in 'pars' is being ignores ad the 'het.source' component of 'ss.opts' is 'FALSE'.")
                pars$sigma.b0.ss <- 0
            }
        }
        if (is.null(ss.link)){
            ss.link <- "identity"
        }
        ## Setting b2.ss to 0 if model is not directional.
        if (!directional){
            pars$b2.ss <- 0
        }
        ## Setting sigma.b0.ss if model does not have heterogeneity in source strengths.
        if (!het.source){
            pars$sigma.b0.ss <- 0
        }
    }
    ## Working out required parameters.
    suppar.names <- c("kappa", "alpha", "sigma.toa")[sim.types[c("bearing", "dist", "toa")]]
    if (sim.ss){
        if (ss.link == "identity"){
            detfn <- "ss"
        } else if (ss.link == "log"){
            detfn <- "log.ss"
        } else if (ss.link == "spherical"){
            detfn <- "spherical.ss"
        } else {
            stop("The argument 'ss.link' must be either \"identity\" or \"log\"")
        }
    }
    detpar.names <- switch(detfn,
                           hn = c("g0", "sigma"),
                           hr = c("g0", "sigma", "z"),
                           th = c("shape", "scale"),
                           lth = c("shape.1", "shape.2", "scale"),
                           ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"),
                           log.ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"),
                           spherical.ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"))
    par.names <- c("D", detpar.names, suppar.names)
    if (!identical(sort(par.names), sort(names(pars)))){
        msg <- paste("The following must be named components of the list 'pars': ",
                     paste(par.names, collapse = ", "), ".", sep = "")
        stop(msg)
    }
    ## Grabbing detection function parameters.
    detpars <- pars[detpar.names]
    out <- vector(mode = "list", length = n.sessions)
    for (s in 1:n.sessions){
        mask <- full.mask[[s]]
        traps <- full.traps[[s]]
        ## Specifies the area in which animal locations can be generated.
        core <- data.frame(x = range(mask[, 1]), y = range(mask[, 2]))
        ## Simulating population.
        if (is.null(cue.rates)){
            if (!popn.provided){
                if (sim.ihd){
                    popn <- as.matrix(sim.popn(D = ihd.surf[[s]], core = mask, buffer = 0,
                                               model2D = "IHP"))
                } else {
                    popn <- as.matrix(sim.popn(D = pars$D, core = core, buffer = 0))
                }
            }
            ## Indicates which individual is being detected.
            individual <- 1:nrow(popn)
        } else {
            D <- pars$D
            if (!first.only){
                ## This is super messy, but it's scaling D from call
                ## density to animal density.
                if (sim.ihd){
                    ihd.surf[[s]] <- ihd.surf[[s]]/mean(cue.rates)
                } else {
                    D <- D/mean(cue.rates)
                }
            }
            if (!popn.provided){
                if (sim.ihd){
                    popn <- as.matrix(sim.popn(D = ihd.surf[[s]], core = mask,
                                               buffer = 0, model2D = "IHP"))
                } else {
                    popn <- as.matrix(sim.popn(D = D, core = core, buffer = 0))
                }
            }
            n.a <- nrow(popn)
            if (freq.dist == "edf"){
                if (length(cue.rates) == 1){
                    freqs <- rep(cue.rates*survey.length, n.a)
                } else {
                    freqs <- sample(cue.rates*survey.length, size = n.a, replace = TRUE)
                }
            } else if (freq.dist == "norm"){
                if (diff(range(cue.rates)) == 0){
                    freqs <- rep(unique(cue.rates)*survey.length, n.a)
                } else {
                    freqs <- rnorm(n.a, mean(cue.rates)*survey.length, sd(cue.rates)*survey.length)
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
            ## Indicates which individual is being detected.
            if (!first.only){
                if (n.a == 0){
                    individual <- numeric(0)
                } else {
                    individual <- rep(1:n.a, times = freqs)
                }
                popn <- popn[individual, , drop=FALSE]
            } else {
                individual <- seq_along(numeric(n.a))
            }
        }
        n.popn <- nrow(popn)
        ## Calculating distances.
        dists <- distances(popn, traps)
        n.traps <- nrow(traps)
        ## Creating empty bincapt if no animals in population.
        if (n.popn == 0){
            captures <- numeric(0)
            bin.capt <- matrix(0, nrow = 0, ncol = n.traps)
            out[[s]] <- list(bincapt = bin.capt)
            if (sim.ss){
                out[[s]]$ss <- bin.capt
            }
        } else {
            ## Calculating detection probabilities and simulating captures.
            if (!sim.ss){
                det.probs <- calc.detfn(dists, detfn, detpars, ss.link)
                if (first.only){
                    ## If only first calls are required, simulate each call separately.
                    full.bin.capt <- matrix(0, nrow = n.a, ncol = n.traps)
                    for (i in 1:n.a){
                        det <- FALSE
                        j <- 1
                        while (!det & j <= freqs[i]){
                            ind.bin.capt <- as.numeric(runif(n.traps) < det.probs[i, ])
                            if (sum(ind.bin.capt) > 0){
                                full.bin.capt[i, ] <- ind.bin.capt
                                det <- TRUE
                            }
                            j <- j + 1
                        }
                    }
                } else {
                    full.bin.capt <- matrix(as.numeric(runif(n.popn*n.traps) < det.probs),
                                            nrow = n.popn, ncol = n.traps)
                }
                captures <- which(apply(full.bin.capt, 1, sum) > 0)
                bin.capt <- full.bin.capt[captures, , drop=FALSE]
                out[[s]] <- list(bincapt = bin.capt)
            } else {
                if (ss.link == "identity"){
                    inv.ss.link <- identity
                } else if (ss.link == "log"){
                    inv.ss.link <- exp
                } else if (ss.link != "spherical"){
                    stop("Argument 'ss.link' must be \"identity\", \"log\", or \"spherical\".")
                }
                pars$cutoff <- cutoff
                detpars$cutoff <- cutoff
                ## Simulating animal directions and calculating orientations
                ## to traps.
                if (pars$b2.ss != 0){
                    if (!is.null(cue.rates) & !first.only){
                        warning("Call directions are being generated independently.")
                    }
                    popn.dirs <- runif(n.popn, 0, 2*pi)
                    popn.bearings <- t(bearings(traps, popn))
                    popn.orientations <- abs(popn.dirs - popn.bearings)
                } else {
                    popn.orientations <- 0
                }
                ## Expected received strength at each microphone for each call.
                if (ss.link %in% c("identity", "log")){
                    ## No idea why something is being divided by 2 in here but I'll assume it's sensible.
                    ss.mean <- inv.ss.link(pars$b0.ss - (pars$b1.ss - pars$b2.ss*(cos(popn.orientations) - 1)/2)*dists)
                } else if (ss.link == "spherical"){
                    ss.mean <- pars$b0.ss - 10*log10(dists^2) - (pars$b1.ss - pars$b2.ss*(cos(popn.orientations) - 1)/2)*(dists - 1)
                }
                ## Random error at each microphone.
                sigma.mat <- matrix(pars$sigma.b0.ss^2, nrow = n.traps, ncol = n.traps)
                diag(sigma.mat) <- diag(sigma.mat) + pars$sigma.ss^2
                if (first.only){
                    if (pars$sigma.b0.ss > 0){
                        stop("Simulation of first call data for situations with heterogeneity in source signal strengths is not yet implemented.")
                        ## Though note that everything is OK for directional calling.
                    }
                    if (inf.calls){
                        log.det.probs <- pnorm(cutoff, ss.mean, pars$sigma.ss, FALSE, TRUE)
                        log.evade.probs <- pnorm(cutoff, ss.mean, pars$sigma.ss, TRUE, TRUE)
                        ## Generating all possible capture histories.
                        n.combins <- 2^n.traps
                        combins <- matrix(NA, nrow = n.combins, ncol = n.traps)
                        for (i in 1:n.traps){
                            combins[, i] <- rep(rep(c(0, 1), each = 2^(n.traps - i)), times = 2^(i - 1))
                        }
                        full.bin.capt <- matrix(0, nrow = n.a, ncol = n.traps)
                        for (i in 1:n.a){
                                        #if (sum(det.probs[i, ]) > 0){
                            log.detprob.mat <- matrix(log.det.probs[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
                            log.evadeprob.mat <- matrix(log.evade.probs[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
                            log.prob.mat <- log.detprob.mat
                            log.prob.mat[combins == 0] <- log.evadeprob.mat[combins == 0]
                            ## Probabilities of each possible capture history.
                            log.d.capt <- apply(log.prob.mat, 1, sum)
                            d.capt <- exp(log.d.capt)
                            ## Selecting a capture history.
                            which.capt <- sample(2:n.combins, size = 1, prob = d.capt[2:n.combins])
                            full.bin.capt[i, ] <- combins[which.capt, ]
                                        #}
                        }
                        full.ss.capt <- full.bin.capt
                        full.ss.capt[full.ss.capt == 1] <- rtruncnorm(sum(full.bin.capt, na.rm = TRUE), a = cutoff,
                                                                      mean = ss.mean[full.bin.capt == 1], sd = pars$sigma.ss)
                    } else {
                        ## If only first calls are required, simulate each call separately.
                        ## Written in C++ as it was way too slow otherwise.
                        full.ss.capt <- sim_ss(ss.mean, pars$sigma.ss, cutoff, freqs)
                    }
                } else {
                    ss.error <- rmvnorm(n.popn, sigma = sigma.mat)
                    ## Filling ss.error for non-hetergeneity models for consistency with old versions.
                    if (pars$sigma.b0.ss == 0){
                        
                        ss.error <- matrix(t(ss.error), nrow = n.popn, ncol = n.traps)
                    }
                    ## Creating SS capture history.
                    full.ss.capt <- ss.mean + ss.error
                }
                captures <- which(apply(full.ss.capt, 1,
                                        function(x, cutoff) any(x > cutoff),
                                        cutoff = cutoff))
                full.bin.capt <- ifelse(full.ss.capt > cutoff, 1, 0)
                ss.capt <- full.ss.capt[captures, , drop=FALSE]
                if (length(captures) == 0){
                    bin.capt <- ss.capt
                } else {
                    bin.capt <- ifelse(ss.capt > cutoff, 1, 0)
                }
                ss.capt[ss.capt < cutoff] <- 0
                out[[s]] <- list(bincapt = bin.capt, ss = ss.capt)
            }
        }
        
        ## Plot to test correct detection simulation.
        if (test.detfn){
            if (!is.null(het.source)){
                if (het.source){
                    warning("Detection function testing for models with heterogeity in source strengths is not yet implemented.")
                    test.detfn <- FALSE
                }
            }
        }
        if (test.detfn & n.popn != 0){
            capt.dists <- dists[full.bin.capt == 1]
            evade.dists <- dists[full.bin.capt == 0]
            all.dists <- c(capt.dists, evade.dists)
            capt.dummy <- c(rep(1, length(capt.dists)),
                            rep(0, length(evade.dists)))
            breaks <- seq(0, max(all.dists), length.out = 1000)
            mids <- breaks[-length(breaks)] + 0.5*diff(breaks)
            breaks[1] <- 0
            split.dummy <- split(capt.dummy,
                                 f = cut(all.dists, breaks = breaks))
            props <- sapply(split.dummy, mean)
            plot(mids, props, type = "l", xlim = c(0, max(all.dists)),
                 ylim = c(0, 1))
            xx <- seq(0, max(all.dists), length.out = 100)
            lines(xx, calc.detfn(xx, detfn, detpars, ss.link), col = "blue")
        }
        ## Total number of detections.
        n.dets <- sum(bin.capt)
        ## Keeping identities of captured individuals.
        capt.ids <- individual[captures]
        ## Locations of captured individuals.
        capt.popn <- popn[captures, ]
        ## IDs of captured individuals.
        ## Capture distances.
        capt.dists <- dists[captures, ]
        ## Simulating additional information.
        if (sim.bearings){
            bearings <- t(bearings(traps, capt.popn))
            bearing.capt <- matrix(0, nrow = nrow(bin.capt),
                                   ncol = ncol(bin.capt))
            bearing.capt[bin.capt == 1] <- (bearings[bin.capt == 1] +
                                            rvm(n.dets, mean = 0, k = pars$kappa)) %% (2*pi)
            out[[s]]$bearing <- bearing.capt
        }
        if (sim.dists){
            dist.capt <- matrix(0, nrow = nrow(bin.capt),
                                ncol = ncol(bin.capt))
            betas <- pars$alpha/capt.dists[bin.capt == 1]
            dist.capt[bin.capt == 1] <- rgamma(n.dets, shape = pars$alpha,
                                               rate = betas)
            out[[s]]$dist <- dist.capt
        }
        if (sim.toas){
            ## Time taken for sound to travel from source to detector.
            toa.capt <- capt.dists/sound.speed*bin.capt
            ## Adding in TOA error.
            toa.capt[bin.capt == 1] <-
                toa.capt[bin.capt == 1] + rnorm(n.dets, sd = pars$sigma.toa)
            out[[s]]$toa <- toa.capt
        }
        if (sim.mrds){
            out[[s]]$mrds <- capt.dists
        }
        if (keep.locs | keep.ids){
            out[[s]] <- list(capt = out[[s]])
            if (keep.locs){
                out[[s]][["capt.locs"]] <- capt.popn
                out[[s]][["popn.locs"]] <- popn
            }
            if (keep.ids){
                out[[s]][["capt.ids"]] <- capt.ids
            }
        }
    }
    if (n.sessions == 1){
        out <- out[[1]]
    }
    out
}

