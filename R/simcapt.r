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
  bincapt <- bincapt[dets, ]
  sim.extra(fit, )
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

#' Simulated supplementary information.
#'
#' Simulates matrices of additional information using parameter estimates
#' from a fitted model returned by \code{admbsecr()}.
#'
#' @param fit a fitted model returnd by \code{\link[admbsecr]{admbsecr}}.
#' @param ... additional arguments
sim.extra <- function(fit, ...){
  UseMethod("sim.extra", fit)
}

#' @S3method sim.extra default
sim.extra.default <- function(fit, ...){
  stop("Method not recognised.")
}

#' @S3method sim.extra simple
sim.extra.simple <- function(fit, bincapt, det.dists, ...){
  bincapt
}

#' @S3method sim.extra ss
sim.extra.ss <- function(fit, bincapt, det.dists, ...){
  bincapt
}

##sim.extra.ang <- function(fit, bincapt, det.dists, ...){
##
##}

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
  coefs <- coef(fit)
  g0 <- ifelse("g0" %in% names(coefs), coefs["g0"], fit$data[["g0"]])
  sigma <- ifelse("sigma" %in% names(coefs), coefs["sigma"], fit$data[["sigma"]])
  g0*exp(-(d^2/(2*sigma^2)))
}

#' @rdname detfn
#' @method detfn hr
#' @S3method detfn hr
detfn.hr <- function(fit, d, ...){
  coefs <- coef(fit)
  g0 <- ifelse("g0" %in% names(coefs), coefs["g0"], fit$data[["g0"]])
  sigma <- ifelse("sigma" %in% names(coefs), coefs["sigma"], fit$data[["sigma"]])
  z <- ifelse("z" %in% names(coefs), coefs["z"], fit$data[["z"]])
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
  coefs <- coef(fit)
  invlink <- c(log = exp, identity = identity)[[link]]
  ssb0 <- ifelse("ssb0" %in% names(coefs), coefs["ssb0"], fit$data[["ssb0"]])
  ssb1 <- ifelse("ssb1" %in% names(coefs), coefs["ssb1"], fit$data[["ssb1"]])
  sigmass <- ifelse("sigmass" %in% names(coefs), coefs["sigmass"], fit$data[["sigmass"]])
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
  coefs <- coef(fit)
  shape <- ifelse("shape" %in% names(coefs), coefs["shape"], fit$data[["shape"]])
  scale <- ifelse("scale" %in% names(coefs), coefs["scale"], fit$data[["scale"]])
  0.5 - 0.5*erf(d/scale - shape)
}
