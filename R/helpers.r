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
getpar <- function(fit, pars){
    allpar.names <- c("D", fit$detpars, fit$suppars)
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
#' @export
sim.capt <- function(fit = NULL, traps = NULL, mask = NULL,
                     infotypes = character(0), detfn = "hn",
                     pars = NULL, ss.link = "identity",
                     cutoff = NULL, sound.speed = 330){
    ## Grabbing values from fit if required.
    ## Specifies the area in which animal locations can be generated.
    core <- data.frame(x = range(mask[, 1]), y = range(mask[, 2]))
    ## Simulating population.
    popn <- sim.popn(D = getpar(fit, "D"), core = core, buffer = 0)
    
}

calc.detfn <- function(d, detfn, pars){
    switch()
}

calc.hn <- function(d, pars){
    if (!identical(sort(names(pars)), c("g0", "sigma"))){
        stop("Components of object 'pars' should have names 'g0' and 'sigma'.")
    }
    g0 <- pars$g0
    sigma <- pars$sigma
    g0*exp(-(d^2/(2*sigma^2)))
}

calc.hr <- function(d, pars){
    if (!identical(sort(names(pars)), c("g0", "sigma", "z"))){
        stop("Components of object 'pars' should have names 'g0', 'sigma', and 'z'.")
    }
    g0 <- pars$g0
    sigma <- pars$sigma
    z <- pars$z
    g0*(1 - exp(-((d/sigma)^-z)))
}

calc.th <- function(d, pars){
    if (!identical(sort(names(pars)), c("scale", "shape"))){
        stop("Components of object 'pars' should have names 'scale' and 'shape'.")
    }
    scale <- pars$scale
    shape <- pars$shape
    0.5 - 0.5*erf(d/scale - shape)
}

calc.logth <- function(d, pars){
    if (!identical(sort(names(pars)), c("scale", "shape1", "shape2"))){
        stop("Components of object 'pars' should have names 'scale' and 'shape1', and 'shape2'.")
    }
    scale <- pars$scale
    shape1 <- pars$shape1
    shape2 <- pars$shape2
    0.5 - 0.5*erf(shape1 - exp(shape2 + scale*d))
}

## Error function
erf <- function(x) 2*pnorm(x*sqrt(2)) - 1
