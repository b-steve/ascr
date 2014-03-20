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
#' @param clicks a data frame containing (at least): (i) \code{$tim},
#' the precise time of arrival of the received sound, and (ii)
#' \code{$trap} the trap at which the sound was recorded.
#' @param dt a \code{K} by \code{K} matrix (where \code{K} is the
#' number of traps) containing the time taken for sound to travel
#' between each pair of traps.
#' @return A data frame. Specifically, the \code{clicks} dataframe,
#' now with a new variable, \code{ID}.
#' @author David Borchers
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

#' Extract estimated or fixed parameter values from an \code{admbsecr} model fit
#'
#' Extracts, estimated, derived, and fitted parameters from a model
#' fitted using \link{admbsecr}.
#'
#' This is a similar function to \link{coef.admbsecr}, however the
#' latter does not allow for extraction of parameters that have been
#' fixed using \code{admbsecr}'s \code{fix} argument.
#'
#' @param pars A character vector containing names of parameter values
#' to be extracted. Alternatively, the character string \code{"all"}
#' will extract all parameters, fixed or otherwise. See the 'Details'
#' section for the \link{admbsecr} function's documentation for
#' information on the parameters that are fitted.
#' @param cutoff Logical, if \code{TRUE}, the cutoff value for an
#' acoustic fit is included.
#' @param as.list Logical, if \code{TRUE}, each parameter value is
#' returned as the component of a list, where component names give the
#' parameter names. In this case, the returned object is suitable as
#' the argument of a variety of functions (e.g., \code{sv} or
#' \code{fix} for the \link{admbsecr} function, or \code{pars} for the
#' \link{sim.capt} function). If \code{FALSE}, parameter values are
#' returned as part of a named vector.
#' @inheritParams locations
#'
#' @return See above information about the argument \code{as.list}. If
#' \code{as.list} is \code{TRUE}, then a list is returned. If
#' \code{as.list} is \code{FALSE}, then a named character vector is
#' returned.
#'
#' @examples
#' get.par(fit = simple.hn.fit, pars = "all")
#' get.par(fit = bearing.hn.fit, pars = c("D", "kappa", "esa"), as.list = TRUE)
#'
#' @export
get.par <- function(fit, pars = "all", cutoff = FALSE, as.list = FALSE){
    allpar.names <- c("D", fit$detpars, fit$suppars, "esa")
    if (length(pars) == 1){
        if (pars == "all"){
            pars <- allpar.names
        } else if (pars == "fitted"){
            pars <- allpar.names[allpar.names != "esa"]
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
    phases$esa <- 0
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
        out <- c(out, fit$cutoff)
        names(out) <- c(pars, "cutoff")
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

#' Extracting mask point locations.
#'
#' Extracts the mask used in an admbsecr fit.
#'
#' @inheritParams locations
#'
#' @return A mask object.
#'
#' @export
get.mask <- function(fit){
    fit$args$mask
}

#' Extracting trap locations.
#'
#' Extracts the trap locations used in an admbsecr fit.
#'
#' @inheritParams locations
#'
#' @return A traps object.
#'
#' @export
get.traps <- function(fit){
    fit$args$traps
}

## Error function.
erf <- function(x){
    2*pnorm(x*sqrt(2)) - 1
}

## Capture probability density surface from a model fit.
p.dot <- function(fit = NULL, points = get.mask(fit), traps = NULL, detfn = NULL,
                  pars = NULL){
    if (!is.null(fit)){
        traps <- get.traps(fit)
        detfn <- fit$args$detfn
        pars <- get.par(fit, fit$detpars, cutoff = fit$fit.types["ss"], as.list = TRUE)
    }
    dists <- distances(traps, points)
    probs <- calc.detfn(dists, detfn, pars)
    aaply(probs, 2, function(x) 1 - prod(1 - x))
}

