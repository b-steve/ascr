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

#' Assigning ID numbers to sounds
#'
#' Identifies recaptures and assigns ID numbers to sounds recorded for
#' an SECR model.
#'
#' Detected sounds are assumed to come from the same animal if times
#' of arrival at different microphones are closer together than the
#' time it would take for sound to travel between these microphones.
#'
#' @param mics a matrix containing the coordinates of trap locations.
#' @param dets a data frame containing (at least): (i) \code{$toa},
#'     the precise time of arrival of the received sound, and (ii)
#'     \code{$trap} the trap at which the sound was recorded.
#' @param sound.speed the speed of sound in metres per second.
#' @return A data frame. Specifically, the \code{dets} dataframe, now
#'     with a new variable, \code{ID}.
#' @author David Borchers
#'
#' @export
make.acoustic.captures <- function(mics, dets, sound.speed){
    mics <- as.matrix(mics)
    dists <- distances(mics, mics)
    dt <- dists/sound.speed
    K <- dim(mics)[1]
    captures <- dets
    ct <- rep(-Inf, K)
    ID <- 1
    ct[dets$trap[1]] <- dets$toa[1]
    new <- FALSE
    ndets <- length(dets$toa)
    for (i in 2:ndets){
        if (ct[dets$trap[i]] > -Inf){
            nd <- length(which(ct > -Inf))
            captures$ID[(i - nd):(i - 1)] <- ID
            ct <- rep(-Inf, K)
            ct[dets$trap[i]] <- dets$toa[i]
            ID <- ID + 1
            if(i == ndets) captures$ID[i] <- ID
        }
        else {
            ct[dets$trap[i]] <- dets$toa[i]
            ctset <- which(ct > -Inf)
            dts <- dt[ctset, dets$trap[i]]
            cts <- -(ct[ctset] - dets$toa[i])
            if (any((cts - dts) > 0)) new <- TRUE
            if (new) {
                nd <- length(which(ct > -Inf)) - 1
                captures$ID[(i - nd):(i - 1)] <- ID
                ct <- rep(-Inf, K)
                ct[dets$trap[i]] <- dets$toa[i]
                ID <- ID + 1
                new <- FALSE
                if (i == ndets) captures$ID[i] <- ID
            } else if(i == ndets){
                nd <- length(which(ct > -Inf))
                captures$ID[(i - nd + 1):i] <- ID
            }
        }
    }
    captures
}

allocate.calls <- function(mics, dets, sound.speed){
    mics <- as.matrix(mics)
    trap.dists <- distances(mics, mics)
    n.dets <- nrow(dets)
    ## Allocating pairwise plausibility of common cue sources.
    dist.mat <- detection_dists(trap.dists, dets$trap)
    timediff.mat <- detection_timediffs(dets$toa, dets$trap)
    maxtime.mat <- dist.mat/sound.speed
    match.mat <- timediff.mat <= maxtime.mat
    ## Finding blocks of multiple cues with possible common sources.
    incomplete.blocks <- find_incomplete_blocks(match.mat)
    n.blocks <- max(incomplete.blocks)
    complete.block <- logical(n.blocks)
    final.mat <- matrix(FALSE, nrow = n.dets, ncol = n.dets)
    ## Allocating possible common cues to sources.
    reqss.mat <- dist.mat/timediff.mat
    for (i in 1:max(incomplete.blocks)){
        ## Grabbing a block.
        block <- match.mat[incomplete.blocks == i, incomplete.blocks == i]
        reqss <- reqss.mat[incomplete.blocks == i, incomplete.blocks == i]
        ## Working out if there is any possible ambiguity.
        is.complete <- all(block)
        ## If ambiguity, resolve it.
        if (!is.complete){
            block <- blockify(block, reqss)
        }
        final.mat[incomplete.blocks == i, incomplete.blocks == i] <- block
    }
    find_incomplete_blocks(final.mat)
}

## Adapted from R2admb.
read.ascr <- function(fn, verbose = FALSE, checkterm = TRUE){
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

#' Extract estimated or fixed parameter values from an \code{ascr} model fit
#'
#' Extracts estimated, derived, and fitted parameters from a model
#' fitted using \link{fit.ascr}.
#'
#' This is a similar function to \link{coef.ascr}, however
#' \code{get.par} also allows for extraction of parameters that have
#' been fixed using \code{fit.ascr}'s \code{fix} argument.
#'
#' @param pars A character vector containing names of parameter values
#'     to be extracted. Alternatively, the character string
#'     \code{"all"} will extract all parameters, fixed or otherwise,
#'     and the character string \code{"fitted"} extracts only fitted
#'     parameters (i.e., not the effective survey area). See the
#'     'Details' section for the \link{fit.ascr} function's
#'     documentation for information on the parameters that are
#'     fitted.
#' @param cutoff Logical, if \code{TRUE}, the cutoff value for an
#'     acoustic fit is included.
#' @param as.list Logical, if \code{TRUE}, each parameter value is
#'     returned as the component of a list, where component names give
#'     the parameter names. In this case, the returned object is
#'     suitable as the argument of a variety of functions (e.g.,
#'     \code{sv} or \code{fix} for the \link{fit.ascr} function, or
#'     \code{pars} for the \link{sim.capt} function). If \code{FALSE},
#'     parameter values are returned as part of a named vector.
#' @inheritParams locations
#'
#' @return See above information about the argument \code{as.list}. If
#'     \code{as.list} is \code{TRUE}, then a list is returned. If
#'     \code{as.list} is \code{FALSE}, then a named character vector
#'     is returned.
#'
#' @examples
#' get.par(fit = example$fits$simple.hn, pars = "all")
#' get.par(fit = example$fits$bearing.hn, pars = c("D", "kappa", "esa"), as.list = TRUE)
#'
#' @export
get.par <- function(fit, pars = "all", cutoff = FALSE, as.list = FALSE){
    esa.names <- paste("esa", 1:fit$n.sessions, sep = ".")
    pars.D <- any((pars == "D" | pars == "all" | pars == "fitted") & !fit$fit.ihd)
    if (pars.D){
        pars[pars == "D"] <- "D.(Intercept)"
    }
    allpar.names <- c(fit$D.betapars, fit$detpars, fit$suppars, esa.names)
    if (length(pars) == 1){
        if (pars == "all"){
            pars <- allpar.names
        } else if (pars == "fitted"){
            pars <- allpar.names[substr(allpar.names, 1, 3) != "esa"]
        }
    }
    if (any(pars == "esa")){
        pars <- pars[-which(pars == "esa")]
        pars <- c(pars, esa.names)
    }
    ## Error checking.
    legal.names <- pars %in% c("D"[pars.D], allpar.names)
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
    if (!fit$fit.types["ss"] & cutoff){
        warning("The cutoff is not being provided as 'fit' does not use signal strength information.")
        cutoff <- FALSE
    }
    out <- numeric(length(pars))
    names(out) <- pars
    det.index <- which(fit$detpars %in% pars)
    supp.index <- which(fit$suppars %in% pars)
    ## Logical vector indicating parameters that weren't estimated.
    phases <- fit$phases
    phases[esa.names] <- 0
    fixed.pars <- phases[pars] == -1
    ## Putting in fixed parameter values.
    if (sum(fixed.pars) > 0){
        out[fixed.pars] <- c(fit$args$sv[pars[fixed.pars]], recursive = TRUE)
    }
    ## Working out parameter groups for parameters in 'pars'.
    det.index <- which(pars %in% fit$detpars)
    supp.index <- which(pars %in% fit$suppars)
    admb.pars <- pars
    ## Putting in estimated parameter values.
    out[!fixed.pars] <- fit$coefficients[admb.pars[!fixed.pars]]
    ## Adding the cutoff if necessary.
    if (cutoff){
        out <- c(out, fit$args$ss.opts$cutoff)
        names(out) <- c(pars, "cutoff")
    }
    if (pars.D){
        names(out)[names(out) == "D.(Intercept)"] <- "D"
        out["D"] <- exp(out["D"])
    }
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

#' Extracting the capture histories
#'
#' Extracts the capture histories from an ascr fit.
#'
#' @inheritParams locations
#' @param session The session from which to extract the capture
#'     histories.
#'
#' @return A capture histories object.
#'
#' @export
get.capt <- function(fit, session = NULL){
    if (fit$n.sessions == 1){
        session <- 1
    }
    if (is.null(session) | !is.list(fit$args$capt)){
        out <- fit$args$capt
    } else {
        out <- fit$args$capt[[session]]
    }
    out
}


#' Extracting mask point locations
#'
#' Extracts the mask used in an ascr fit.
#'
#' @inheritParams locations
#' @param session The session from which to extract the mask.
#'
#' @return A mask object.
#'
#' @export
get.mask <- function(fit, session = NULL){
    if (fit$n.sessions == 1){
        session <- 1
    }
    if (is.null(session) | !is.list(fit$args$mask)){
        out <- fit$args$mask
    } else {
        out <- fit$args$mask[[session]]
    }
    out
}

#' Extracting trap locations
#'
#' Extracts the trap locations used in an ascr fit.
#'
#' @inheritParams locations
#' @param session The session from which to extract the trap
#'     locations.
#' 
#' @return A traps object.
#' 
#' @export
get.traps <- function(fit, session = NULL){
    if (fit$n.sessions == 1){
        session <- 1
    }
    if (is.null(session) | !is.list(fit$args$traps)){
        out <- fit$args$traps
    } else {
        out <- fit$args$traps[[session]]
    }
    out
}

## Error function.
erf <- function(x){
    2*pnorm(x*sqrt(2)) - 1
}

#' Calculating detection probabilities.
#'
#' Calculates the probability of detection by at least one detector
#' for specific locations in the survey area.
#'
#' @param fit A fitted model from \link{fit.ascr}.
#' @param session For multisession models, the session from which the
#'     trap locations should be taken.
#' @param esa Logical, if \code{TRUE} the effective sampling area is
#'     returned instead of capture probabilities.
#' @param points A matrix with two columns. Each row provides
#'     Cartesian coordinates for the location of a point at which a
#'     capture probability should be returned.
#' @param traps A matrix with two columns. Each row provides Cartesian
#'     coordinates for the location of a trap (or detector). Ignored
#'     if \code{fit} is not \code{NULL}.
#' @param detfn A character string specifying the detection function
#'     to be used. One of "hn" (halfnormal), "hr" (hazard rate), "th"
#'     (threshold), "lth" (log-link threshold), or "ss" (signal
#'     strength). Ignored if \code{fit} is not \code{NULL}.
#' @param ss.link A character string, either \code{"identity"},
#'     \code{"log"}, or \code{"spherical"}, which specifies the
#'     relationship between the expected received signal strength and
#'     distance from the microphone. See the documentation for
#'     \link{fit.ascr} for further details. Ignored if \code{fit} is
#'     not \code{NULL}.
#' @param pars A named list. Component names are parameter names, and
#'     each component is a value for the associated parameter. Ignored
#'     if \code{fit} is not \code{NULL}.
#' @param n.quadpoints An integer, giving the number of quadrature
#'     points used for numerical integration over the possible call
#'     directions.
#'
#' @return A vector containing detection probabilities for each
#'     location in \code{points}.
#'
#' @export
p.dot <- function(fit = NULL, session = 1, esa = FALSE, points = get.mask(fit, session),
                  traps = NULL, detfn = NULL, ss.link = NULL, pars = NULL, n.quadpoints = 8){
    if (!is.null(fit)){
        if (is.null(traps)){
            traps <- get.traps(fit, session)
        }
        detfn <- fit$args$detfn
        pars <- get.par(fit, fit$detpars, cutoff = fit$fit.types["ss"], as.list = TRUE)
        ss.link <- fit$args$ss.opts$ss.link
        re.detfn <- fit$re.detfn
    } else {
        re.detfn <- FALSE
        if (detfn == "ss"){
            if (pars$b2.ss != 0 | pars$sigma.b0.ss != 0){
                re.detfn <- TRUE
            }
        }
    }
    dists <- distances(traps, points)
    ## Calculating probabilities of detection when random effects are
    ## in detection function. Detections at traps no longer
    ## independent.
    if (re.detfn){
        if (!is.null(pars$sigma.b0.ss)){
            if (pars$sigma.b0.ss != 0){
                stop("Function p.dot() has not yet been implemented for models with heterogeneous source strengths.")
            }
        }
        n.traps <- nrow(traps)
        n.points <- nrow(points)
        dirs <- (0:(n.quadpoints - 1))*2*pi/n.quadpoints
        probs <- numeric(n.points)
        ## Integrating over all possible directions.
        ## TODO: Write all this in C++.
        for (i in 1:n.quadpoints){
            dir <- dirs[i]
            bearings <- bearings(traps, points)
            orientations <- abs(dir - bearings)
            for (j in 1:n.points){
                ## Probabilities of detection given orientation.
                o.prob <- numeric(n.traps)
                for (k in 1:n.traps){
                    o.prob[k] <- calc.detfn(dists[k, j], detfn, pars, ss.link,
                                            orientations[k, j])
                }
                probs[j] <- probs[j] + (1/n.quadpoints)*(1 - prod(1 - o.prob))
            }
        }
        out <- probs
    } else {
        probs <- calc.detfn(dists, detfn, pars, ss.link)
        out <- aaply(probs, 2, function(x) 1 - prod(1 - x))
    }
    if (esa){
        A <- attr(points, "area")
        if (is.null(A)){
            stop("The argument 'points' must be a mask object if ESA is to be calculated.")
        }
        out <- A*sum(out)
    }
    out
}

## Link functions for fit.ascr() function.
log.link <- function(x){
    x <- pmax(x, .Machine$double.eps)
    log(x)
}

logit.link <- function(x){
    x <- pmax(x, .Machine$double.eps)
    x <- pmin(x, 1 - .Machine$double.neg.eps)
    log(x/(1 - x))
}

scaled.logit.link <- function(x){
    (x/10)^3
}

inv.logit <- function(x){
    exp(x)/(exp(x) + 1)
}

inv.scaled.logit.link <- function(x){
    x <- inv.logit(x)
    2*1e8*x - 1e8
}

scaled.log.link <- function(x){
    x <- x + 1e8
    log(x)
}

#' Extracting Monte Carlo error
#'
#' Extracts calculated Monte Carlo errors from a bootstrap procedure
#' carried out by \link{boot.ascr}.
#'
#' @param estimate A character string, either \code{"bias"} or
#'     \code{"se"}, which determines whether Monte Carlo errors for
#'     bias estimates or standard errors are reported.
#' @inheritParams locations
#'
#' @seealso \link{boot.ascr} for the bootstrap procedure.
#' @seealso \link{stdEr.ascr.boot} for standard errors.
#' @seealso \link{get.bias} for estimated biases.
#'
#' @export
get.mce <- function(fit, estimate){
    if (estimate == "bias"){
        out <- fit$boot$bias.mce
    } else if (estimate == "se"){
        out <- fit$boot$se.mce
    } else {
        stop("The argument 'estimate' must be either \"bias\" or \"se\"")
    }
    out
}

#' Extracting estimated biases
#'
#' Extracts bias in parameter estimates, calculated using the
#' bootstrap procedure carried out by \link{boot.ascr}.
#'
#' @inheritParams locations
#' @inheritParams coef.ascr
#' @inheritParams stdEr.ascr.boot
#'
#' @seealso \link{boot.ascr} for the bootstrap procedure.
#' @seealso \link{get.mce} for Monte Carlo error the biases are
#'     subject to.
#'
#' @export
get.bias <- function(fit, pars = "fitted", mce = FALSE){
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    if (any(pars == "esa")){
        pars <- pars[-which(pars == "esa")]
        pars <- c(pars, paste("esa", 1:fit$n.sessions, sep = "."))
    }
    par.names <- names(fit$coefficients)
    if (!all(pars %in% c("fitted", "derived", "linked", par.names))){
        stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
    }
    mces <- get.mce(fit, estimate = "bias")
    if (any(c("fitted", "derived", "linked") %in% pars)){
        which.linked <- grep("_link", par.names)
        linked <- fit$boot$bias[which.linked]
        which.derived <- which(substr(par.names, 1, 3) == "esa" | par.names == "Da")
        derived <- fit$boot$bias[which.derived]
        fitted <- fit$boot$bias[-c(which.linked, which.derived)]
        out <- mget(pars)
        names(out) <- NULL
        out <- c(out, recursive = TRUE)
        if (!fit$fit.ihd){
            out <- out[c("D", names(out)[!(names(out) %in% c("D.(Intercept)", "D"))])]
        }
    } else {
        out <- fit$boot$bias[pars]
    }
    if (mce){
        out.vec <- out
        out <- cbind(out.vec, mces[names(out)])
        rownames(out) <- names(out.vec)
        colnames(out) <- c("Bias", "MCE")
    }
    out
}

## Helper function for passing lists as 3D arrays to ADMB.
list.to.vector <- function(x){
    c(lapply(x, function(m) c(t(m))), recursive = TRUE)
}

## Calculates the `effective listening area', assuming perfect
## detection within some fixed radius of the traps.
calc.ela <- function(traps, radius, mask = NULL, ...){
    if (is.null(mask)){
        mask <- create.mask(traps, buffer = 1.1*radius, ...)
    }
    a <- attr(mask, "area")
    in.area <- apply(distances(mask, as.matrix(traps)), 1, function(x) min(x) < radius)
    a*sum(in.area)
}
