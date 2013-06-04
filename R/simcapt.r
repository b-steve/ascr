## Simulates capt matrix suitable for admbsecr(), from either a model
## fit (using estimated parameters) or from provided values.

sim.capt <- function(fit){
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
  buffer <- 0
  range.x <- range(mask[, 1])
  range.y <- range(mask[, 2])
  core <- data.frame(x = range.x, y = range.y)
  popn <- sim.popn(D = allcoefs["D"], core = core, buffer = buffer)
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

#' @S3method sim.extra ss
sim.extra.ss <- function(fit, bincapt, det.dists, det.locs, ...){
  dim <- dim(bincapt)
  dim(bincapt) <- c(dim[1], 1, dim[2])
  bincapt
}

#' @S3method sim.extra ang
sim.extra.ang <- function(fit, bincapt, det.dists, det.locs, ...){
  traps <- fit$traps
  kappa <- getpar(fit, "kappa")
  dim <- dim(bincapt)
  ndets <- sum(bincapt)
  angcapt <- angles(det.locs, traps)*bincapt
  angcapt[bincapt == 1] <- (angcapt[bincapt == 1] +
      rvm(ndets, mean = 1, k = kappa)) %% (2*pi)
  angcapt[angcapt == 0 & bincapt == 1] <- 2*pi
  dim(angcapt) <- c(dim[1], 1, dim[2])
  angcapt
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
  angcapt <- angles(det.locs, traps)*bincapt
  angcapt[bincapt == 1] <- (angcapt[bincapt == 1] +
      rvm(ndets, mean = 1, k = kappa)) %% (2*pi)
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
