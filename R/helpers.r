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
#' @param fit A fitted model from \link[admbsecr]{admbsecr}.
#' @param type A character string. If \code{"fixed"} then only
#' standard errors for estimated parameters are presented. If
#' \code{"all"} then standard errors are also presented for derived
#' paramters (e.g., the effective survey area).
#' @export
stdEr <- function(fit, type = "fixed"){
    out <- fit$se
    if (type == "fixed"){
        out <- out[-length(out)]
    }
    out
}

## Return fixed or estimated parameter values from a model fit.
getpar <- function(fit, pars, cutoff = FALSE, as.list = FALSE){
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
    if (!fit$fit.types["ss"] & cutoff){
        warning("The cutoff is not being provided as 'fit' does not use signal strength information.")
        cutoff <- FALSE
    }
    out <- numeric(length(pars))
    names(out) <- pars
    det.index <- which(fit$detpars %in% pars)
    supp.index <- which(fit$suppars %in% pars)
    ## Logical vector indicating parameters that weren't estimated.
    fixed.pars <- fit$phases[pars] == -1
    ## Putting in fixed parameter values.
    if (sum(fixed.pars) > 0){
        out[fixed.pars] <- c(fit$sv[pars[fixed.pars]], recursive = TRUE)
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

## Error function.
erf <- function(x){
    2*pnorm(x*sqrt(2)) - 1
}

## Capture probability density surface from a model fit.
p.dot <- function(fit = NULL, points = fit$mask, traps = NULL, detfn = NULL,
                  pars = NULL){
    if (!is.null(fit)){
        traps <- fit$traps
        detfn <- fit$detfn
        pars <- getpar(fit, fit$detpars, cutoff = fit$fit.types["ss"], as.list = TRUE)
    }
    dists <- distances(traps, points)
    probs <- calc.detfn(dists, detfn, pars)
    aaply(probs, 2, function(x) 1 - prod(1 - x))
}
