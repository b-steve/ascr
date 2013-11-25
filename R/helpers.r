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
#' Identifies recaptures and assigns ID numbers to sounds recorded for an SECR model.
#'
#' Detected sounds are assumed to come from the same animal if times of arrival at
#' different microphones are closer together than the time it would take for sound to
#' travel between these microphones.
#'
#' @param mics a matrix containing the coordinates of trap locations.
#' @param clicks a data frame containing (at least): (i) \code{$tim$}, the precise
#' time of arrival of the received sound, and (ii) \code{$trap} the trap at which
#' the sound was recorded.
#' @param dt a \code{K} by \code{K} matrix (where \code{K} is the number of traps)
#' containing the time taken for sound to travel between each pair of traps.
#' @return A data frame. Specifically, the \code{clicks} dataframe, now with a new
#' variable, \code{ID}.
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
#' @param type a character string, either \code{"fixed"}, or
#' \code{"all"}. If \code{"fixed"} only model parameter standard errors
#' are shown, otherwise a standard error (calculated using the delta
#' method) is also provided for the effective sampling area.
#' @export
stdEr <- function(fit, type = "fixed"){
  out <- fit$se
  if (type == "fixed"){
    out <- out[names(out) != "esa"]
  }
}
