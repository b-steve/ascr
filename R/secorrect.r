#' SECR standard error correction
#'
#' Corrects standard errors from an \code{admbsecr} fit.
#'
#' @param fit an \code{admbsecr} model fit.
#' @param size number of bootstrap samples.
#' @param calls a vector containing the numbers of calls emitted from
#' recorded individuals. If \code{NULL}, this defaults to what
#' is provided by \code{fit}.
#' @export
se.correct <- function(fit, size, calls = NULL){
  coefs <- coef(fit, type = "all")
  start <- coef(fit)
  npars <- length(coefs)
  res <- matrix(0, nrow = size, ncol = npars + 1)
  traps <- fit[["traps"]]
  mask <- fit[["mask"]]
  bounds <- fit[["bounds"]]
  fix <- fit[["fix"]]
  cutoff <- fit$data[["c"]]
  if (is.null(calls)){
    cpi <- fit$data[["cpi"]]
  } else {
    cpi <- calls
  }
  boot.calls <- length(cpi) > 1
  sound.speed <- fit[["sound.speed"]]
  method <- fit[["method"]]
  detfn <- fit[["detfn"]]
  scalefactors <- fit[["scalefactors"]]
  memory <- fit[["memory"]]
  colnames(res) <- c(names(coefs), "maxgrad")
  for (i in 1:size){
    if (boot.calls){
      cpi.boot <- sample(cpi, replace = TRUE)
    } else {
      cpi.boot <- NULL
    }
    capt <- sim.capt(fit = fit, calls = round(cpi))
    bootfit <- try.admbsecr(sv = start, capt = capt, traps = traps, mask = mask,
                            bounds = bounds, fix = fix, cutoff = cutoff,
                            cpi = cpi.boot, sound.speed = sound.speed,
                            method = method, detfn = detfn, memory = memory,
                            scalefactors = scalefactors)
    if (class(bootfit)[1] == "try-error"){
      res[i, ] <- NA
      warning(paste("Failed convergence on iteration", i, sep = " "))
    } else {
      res[i, ] <- c(coef(bootfit, type = "all"), bootfit$maxgrad)
    }
  }
  conv <- res[, "maxgrad"] > -1
  res <- res[conv, -which(colnames(res) == "maxgrad")]
  bias <- apply(res, 2, mean, na.rm = TRUE) - coefs
  coefficients.corrected <- coefs - bias
  se.corrected <- apply(res, 2, sd, na.rm = TRUE)
  out <- list(boots = res, coefficients.corrected = coefficients.corrected,
              se.corrected = se.corrected)
  fit$se.correct <- out
  class(fit) <- c("secorrect", class(fit))
  fit
}

## Modified from R2admb.
#' @S3method summary secorrect
summary.secorrect <- function(object, ...){
  parnames <- names(object$se.correct$coefficients.corrected)
  coef.p <- object$se.correct$coefficients.corrected[parnames != "D"]
  s.err <- object$se.correct$se.corrected[parnames != "D"]
  tvalue <- coef.p/s.err
  dn <- c("Estimate", "Std. Error")
  pvalue <- 2 * pnorm(-abs(tvalue))
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                "z value", "Pr(>|z|)"))
  ans <- c(list(coefficients=coef.table),
           object[c("loglik","fn","npar")])
  class(ans) <- "summary.admb"
  ans
}

#' @S3method coef secorrect
coef.secorrect <- function(object, ...){
  object$se.correct$coefficients.corrected
}
