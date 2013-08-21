## Simulates capt matrix suitable for admbsecr(), from either a model
## fit (using estimated parameters) or from provided values.

#' Simulating SECR data
#'
#' Simulates SECR data capture histories in the correct format for use with
#' \code{\link[admbsecr]{admbsecr}}.
#'
#' If \code{fit} is provided then no other arguments are required. Otherwise, at least
#' \code{traps}, \code{mask}, \code{pars}, \code{method}, and \code{detfn} are needed.
#'
#' @param fit a fitted \code{admbsecr} model object which provides the method,
#' detection function, and parameter values from which to generate capture
#' histories.
#' @param calls a vector containing call frequencies from monitored individuals.
#' @param traps a matrix containing the coordinates of trap locations. The object
#' returned by \code{\link[secr]{read.traps}} is suitable.
#' @param mask a mask object. The object returned by \code{\link[secr]{make.mask}} is
#' suitable.
#' @param pars a named vector of parameter values.
#' @param method the \code{admbsecr} method which specifies the additional information
#' to be simulated.
#' @param detfn the admbsecr detection function to be used.
#' @param cutoff the signal strength threshold of detection. Required if \code{method}
#' is \code{"ss"} or \code{"sstoa"}.
#' @param sound.speed the speed of sound in metres per second. Used for TOA analysis.
#' @export
sim.capt <- function(fit, calls = 1, traps = NULL, mask = NULL, pars = NULL,
                     method = "simple", detfn = "hn", cutoff = NULL,
                     sound.speed = NULL){
  if (!missing(fit)){
    method <- fit$method
    detfn <- fit$detfn
    traps <- gettraps(fit)
    mask <- getmask(fit)
    allparnames <- fit$parnames
    allcoefs <- getpar(fit, allparnames)
  } else {
    allcoefs <- pars
    fit <- list(traps = traps, mask = mask, coefficients = allcoefs,
                data = list(c = cutoff), method = method, detfn = detfn)
    class(fit) <- c("admbsim", method, detfn)
  }
  ntraps <- nrow(traps)
  buffer <- 0
  range.x <- range(mask[, 1])
  range.y <- range(mask[, 2])
  core <- data.frame(x = range.x, y = range.y)
  popn <- as.matrix(sim.popn(D = getpar(fit, "D")/mean(calls),
                             core = core, buffer = buffer))
  popn.dists <- distances(as.matrix(popn), as.matrix(traps))
  if (length(calls) != 1 | calls[1] != 1){
    n <- nrow(popn)
    freqs <- resample(calls, size = n, replace = TRUE)
    ncalls <- sum(freqs)
    pos <- c(1, cumsum(freqs) + 1)
    call.popn <- matrix(0, nrow = ncalls, ncol = 2)
    rownames(call.popn) <- 1:ncalls
    colnames(call.popn) <- c("x", "y")
    call.popn.dists <- matrix(0, nrow = ncalls, ncol = ntraps)
    rownames(call.popn.dists) <- 1:ncalls
    colnames(call.popn.dists) <- 1:ntraps
    for (i in 1:n){
      for (j in 1:2){
        call.popn[pos[i]:(pos[i + 1] - 1), j] <- popn[i, j]
      }
      for (j in 1:ntraps){
        call.popn.dists[pos[i]:(pos[i + 1] - 1), j] <- popn.dists[i, j]
      }
    }
    popn <- call.popn
    popn.dists <- call.popn.dists
  }
  bincapt <- sim.bincapt(fit, popn.dists)
  dets <- apply(bincapt, 1, function(x) any(x != 0))
  det.dists <- popn.dists[dets, ]
  det.locs <- popn[dets, ]
  bincapt <- bincapt[dets, ]
  sim.extra(fit, bincapt, det.dists, det.locs)
}

## Will produce a binary capture history for all methods other than
## those using a signal strength detection function.
sim.bincapt <- function(fit, popn.dists){
  dim <- dim(popn.dists)
  if (fit$method == "ss" | fit$method == "sstoa"){
    coefs <- coef(fit)
    cutoff <- fit$data[["c"]]
    muss <- calc.detfn(fit, popn.dists, prob = FALSE)
    sigmass <- ifelse("sigmass" %in% names(coefs), coefs["sigmass"],
                      fit$data[["sigmass"]])
    errors <- matrix(rnorm(prod(dim), sd = sigmass),
                     nrow = dim[1], ncol = dim[2])
    ss <- muss + errors
    bincapt <- ifelse(ss < cutoff, 0, ss)
  } else {
    capt.probs <- calc.detfn(fit, popn.dists)
    capt.rnos <- matrix(runif(prod(dim)), nrow = dim[1], ncol = dim[2])
    bincapt <- capt.rnos <= capt.probs
    bincapt <- ifelse(bincapt, 1, 0)
  }
  bincapt
}

sim.extra <- function(fit, ...){
  UseMethod("sim.extra", fit)
}

#' @S3method sim.extra default
sim.extra.default <- function(fit, ...){
  stop("Method not recognised.")
}

#' @S3method sim.extra simple
sim.extra.simple <- function(fit, bincapt, det.dists, det.locs, ...){
  dim <- dim(bincapt)
  dim(bincapt) <- c(dim[1], 1, dim[2])
  bincapt
}

#' @S3method sim.extra toa
sim.extra.toa <- function(fit, bincapt, det.dists, det.locs, ...){
  sigmatoa <- getpar(fit, "sigmatoa")
  dim <- dim(bincapt)
  n <- dim[1]
  ntraps <- dim[2]
  ndets <- sum(bincapt)
  errors <- rnorm(ndets, 0, sigmatoa)
  ## Start times are arbitrarily 100*id.no
  toacapt <- matrix(rep(100*(1:n), ntraps), nrow = n, ncol = ntraps)*bincapt
  ## Need to set up expected arrival times
  travel.time <- det.dists/330
  toacapt[bincapt == 1] <- toacapt[bincapt == 1] + travel.time[bincapt == 1] +
      errors
  dim(toacapt) <- c(dim[1], 1, dim[2])
  toacapt
}

#' @S3method sim.extra ang
sim.extra.ang <- function(fit, bincapt, det.dists, det.locs, ...){
  traps <- gettraps(fit)
  kappa <- getpar(fit, "kappa")
  dim <- dim(bincapt)
  ndets <- sum(bincapt)
  angcapt <- t(angles(traps, det.locs))*bincapt
  angcapt[bincapt == 1] <- (angcapt[bincapt == 1] +
      rvm(ndets, mean = 0, k = kappa)) %% (2*pi)
  angcapt[angcapt == 0 & bincapt == 1] <- 2*pi
  dim(angcapt) <- c(dim[1], 1, dim[2])
  angcapt
}

## Does not do anything as bincapt already returns SS information.
#' @S3method sim.extra ss
sim.extra.ss <- function(fit, bincapt, det.dists, det.locs, ...){
  dim <- dim(bincapt)
  dim(bincapt) <- c(dim[1], 1, dim[2])
  bincapt
}

## Does not do anything as bincapt already returns SS information.
#' @S3method sim.extra sstoa
sim.extra.sstoa <- function(fit, bincapt, det.dists, det.locs, ...){
  sscapt <- bincapt
  ## Add zero to coerce to numeric.
  bincapt <- 0 + bincapt != 0
  sigmatoa <- getpar(fit, "sigmatoa")
  dim <- dim(bincapt)
  n <- dim[1]
  ntraps <- dim[2]
  ndets <- sum(bincapt)
  errors <- rnorm(ndets, 0, sigmatoa)
  ## Start times are arbitrarily 100*id.no
  toacapt <- matrix(rep(100*(1:n), ntraps), nrow = n, ncol = ntraps)*bincapt
  ## Need to set up expected arrival times
  travel.time <- det.dists/330
  toacapt[bincapt == 1] <- toacapt[bincapt == 1] + travel.time[bincapt == 1] +
      errors
  dim(toacapt) <- c(dim[1], 1, dim[2])
  dim(sscapt) <- c(dim[1], 1, dim[2])
  array(c(sscapt, toacapt), dim = c(dim[1], 1, dim[2], 2))
}

#' @S3method sim.extra dist
sim.extra.dist <- function(fit, bincapt, det.dists, det.locs, ...){
  alpha <- getpar(fit, "alpha")
  dim <- dim(bincapt)
  ndets <- sum(bincapt)
  distcapt <- det.dists*bincapt
  betas <- alpha/distcapt[bincapt == 1]
  distcapt[bincapt == 1] <- rgamma(ndets, shape = alpha, rate = betas)
  dim(distcapt) <- c(dim[1], 1, dim[2])
  distcapt
}

#' @S3method sim.extra angdist
sim.extra.angdist <- function(fit, bincapt, det.dists, det.locs, ...){
  traps <- gettraps(fit)
  alpha <- getpar(fit, "alpha")
  kappa <- getpar(fit, "kappa")
  dim <- dim(bincapt)
  ndets <- sum(bincapt)
  angcapt <- t(angles(traps, det.locs))*bincapt
  angcapt[bincapt == 1] <- (angcapt[bincapt == 1] +
      rvm(ndets, mean = 0, k = kappa)) %% (2*pi)
  angcapt[angcapt == 0 & bincapt == 1] <- 2*pi
  dim(angcapt) <- c(dim[1], 1, dim[2])
  distcapt <- det.dists*bincapt
  betas <- alpha/distcapt[bincapt == 1]
  distcapt[bincapt == 1] <- rgamma(ndets, shape = alpha, rate = betas)
  dim(distcapt) <- dim(angcapt)
  jointcapt <- array(0, dim = c(dim(angcapt), 2))
  jointcapt[, , , 1] <- angcapt
  jointcapt[, , , 2] <- distcapt
  jointcapt
}

#' @S3method sim.extra mrds
sim.extra.mrds <- function(fit, bincapt, det.dists, det.locs, ...){
  ## TO DO
  stop("Method not yet implemented.")
}

#' @S3method sim.extra ssmrds
sim.extra.ssmrds <- function(fit, bincapt, det.dists, det.locs, ...){
  ## TO DO
  stop("Method not yet implemented.")
}

#' Estimated detection probability from a fitted model.
#'
#' Calculates the probability of detection for given distances from a
#' fit returned by \code{admbsecr()}.
#'
#' @param fit a fitted model returned by \code{\link[admbsecr]{admbsecr}}.
#' @param ... additional arguments.
#' @rdname calc.detfn
#' @export
calc.detfn <- function(fit, ...){
  UseMethod("calc.detfn", fit)
}

#' @S3method calc.detfn default
calc.detfn.default <- function(fit, ...){
  stop(paste("This function does not work with admbsecr fits with detection function ",
             "\"", fit$calc.detfn, "\"", sep = ""))
}

#' @rdname calc.detfn
#' @param d vector of distances from which probabilities are calculated.
#' @method calc.detfn hn
#' @S3method calc.detfn hn
calc.detfn.hn <- function(fit, d, ...){
  g0 <- getpar(fit, "g0")
  sigma <- getpar(fit, "sigma")
  g0*exp(-(d^2/(2*sigma^2)))
}

#' @rdname calc.detfn
#' @method calc.detfn hr
#' @S3method calc.detfn hr
calc.detfn.hr <- function(fit, d, ...){
  g0 <- getpar(fit, "g0")
  sigma <- getpar(fit, "sigma")
  z <- getpar(fit, "z")
  g0*(1 - exp(-((d/sigma)^-z)))
}

## Returns probabilities if prob = TRUE, otherwise returns expected
## signal strengths.
#' @rdname calc.detfn
#' @param prob logical, if \code{TRUE}, capture probability is returned.
#' If \code{FALSE}, expected signal strengths are returned.
#' @method  calc.detfn ss
#' @S3method calc.detfn ss
calc.detfn.ss <- function(fit, d, prob = TRUE, ...){
  link <- fit$detfn
  cutoff <- fit$data[["c"]]
  invlink <- c(log = exp, identity = identity)[[link]]
  ssb0 <- getpar(fit, "ssb0")
  ssb1 <- getpar(fit, "ssb1")
  sigmass <- getpar(fit, "sigmass")
  muss <- invlink(ssb0 - ssb1*d)
  if (prob){
    out <- 1 - pnorm((cutoff - muss)/sigmass)
  } else {
    out <- muss
  }
  out
}

#' @rdname calc.detfn
#' @method  calc.detfn sstoa
#' @S3method calc.detfn sstoa
calc.detfn.sstoa <- calc.detfn.ss

#' @rdname calc.detfn
#' @method calc.detfn th
#' @S3method calc.detfn th
calc.detfn.th <- function(fit, d, ...){
  shape <- getpar(fit, "shape")
  scale <- getpar(fit, "scale")
  0.5 - 0.5*erf(d/scale - shape)
}

#' @S3method coef admbsim
coef.admbsim <- function(object, ...){
  object$coefficients
}

resample <- function(x, ...){
  x[sample.int(length(x), ...)]
}
