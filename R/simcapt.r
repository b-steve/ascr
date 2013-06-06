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
#' @param calls number of calls to be emitted by each individual.
#' @export
sim.capt <- function(fit, traps = NULL, mask = NULL, pars = NULL, method = "simple",
                     detfn = "hn", cutoff = NULL, sound.speed = NULL, calls = 1){
  if (!missing(fit)){
    method <- fit$method
    detfn <- fit$detfn
    traps <- fit$traps
    mask <- fit$mask
    fitcoefs <- coef(fit)
    allparnames <- fit$parnames
    coefparnames <- names(fitcoefs)
    dataparnames <- allparnames[!allparnames %in% coefparnames]
    fixcoefs <- numeric(length(dataparnames))
    names(fixcoefs) <- dataparnames
    for (i in dataparnames){
      fixcoefs[i] <- fit$data[[i]]
    }
    allcoefs <- c(fitcoefs, fixcoefs)
  } else {
    allcoefs <- pars
    fit <- list(traps = traps, mask = mask, coefficients = allcoefs,
                data = list(c = cutoff), method = method, detfn = detfn)
    class(fit) <- c("admbsim", method, detfn)
  }
  buffer <- 0
  range.x <- range(mask[, 1])
  range.y <- range(mask[, 2])
  core <- data.frame(x = range.x, y = range.y)
  popn <- sim.popn(D = allcoefs["D"]/calls, core = core, buffer = buffer)
  if (calls != 1){
    popn <- 0 ## UNFINISHED
  }
  popn.dists <- distances(as.matrix(popn), as.matrix(traps))
  bincapt <- sim.bincapt(fit, popn.dists)
  dets <- apply(bincapt, 1, function(x) any(x != 0))
  det.dists <- popn.dists[dets, ]
  det.locs <- popn[dets, ]
  bincapt <- bincapt[dets, ]
  sim.extra(fit, bincapt, det.dists, det.locs)
}

sim.bincapt <- function(fit, popn.dists){
  dim <- dim(popn.dists)
  if (fit$method == "ss" | fit$method == "sstoa"){
    coefs <- coef(fit)
    cutoff <- fit$data[["c"]]
    muss <- detfn(fit, popn.dists, prob = FALSE)
    sigmass <- ifelse("sigmass" %in% names(coefs), coefs["sigmass"],
                      fit$data[["sigmass"]])
    errors <- matrix(rnorm(prod(dim), sd = sigmass),
                     nrow = dim[1], ncol = dim[2])
    ss <- muss + errors
    bincapt <- ifelse(ss < cutoff, 0, ss)
  } else {
    capt.probs <- detfn(fit, popn.dists)
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
  ## TO DO
  stop("Method not yet implemented.")
}

#' @S3method sim.extra ang
sim.extra.ang <- function(fit, bincapt, det.dists, det.locs, ...){
  traps <- fit$traps
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

#' @S3method sim.extra ss
sim.extra.ss <- function(fit, bincapt, det.dists, det.locs, ...){
  dim <- dim(bincapt)
  dim(bincapt) <- c(dim[1], 1, dim[2])
  bincapt
}

#' @S3method sim.extra sstoa
sim.extra.sstoa <- function(fit, bincapt, det.dists, det.locs, ...){
  ## TO DO
  stop("Method not yet implemented.")
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
  traps <- fit$traps
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
#' @rdname detfn
#' @export
detfn <- function(fit, ...){
  UseMethod("detfn", fit)
}

#' @S3method detfn default
detfn.default <- function(fit, ...){
  stop(paste("This function does not work with admbsecr fits with detection function ",
             "\"", fit$detfn, "\"", sep = ""))
}

#' @rdname detfn
#' @param d vector of distances from which probabilities are calculated.
#' @method detfn hn
#' @S3method detfn hn
detfn.hn <- function(fit, d, ...){
  g0 <- getpar(fit, "g0")
  sigma <- getpar(fit, "sigma")
  g0*exp(-(d^2/(2*sigma^2)))
}

#' @rdname detfn
#' @method detfn hr
#' @S3method detfn hr
detfn.hr <- function(fit, d, ...){
  g0 <- getpar(fit, "g0")
  sigma <- getpar(fit, "sigma")
  z <- getpar(fit, "z")
  g0*(1 - exp(-((d/sigma)^-z)))
}

## Returns probabilities if prob = TRUE, otherwise returns expected
## signal strengths.
#' @rdname detfn
#' @param prob logical, if \code{TRUE}, capture probability is returned.
#' If \code{FALSE}, expected signal strengths are returned.
#' @method  detfn ss
#' @S3method detfn ss
detfn.ss <- function(fit, d, prob = TRUE, ...){
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

#' @rdname detfn
#' @method  detfn sstoa
#' @S3method detfn sstoa
detfn.sstoa <- detfn.ss

#' @rdname detfn
#' @method detfn th
#' @S3method detfn th
detfn.th <- function(fit, d, ...){
  shape <- getpar(fit, "shape")
  scale <- getpar(fit, "scale")
  0.5 - 0.5*erf(d/scale - shape)
}

#' @S3method coef admbsim
coef.admbsim <- function(object, ...){
  object$coefficients
}
