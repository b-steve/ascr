#' SECR standard error correction
#'
#' Corrects standard errors from an \code{admbsecr} fit.
#'
#' @param fit an \code{admbsecr} model fit.
#' @param calls a vector containing call frequencies from monitored individuals.
#' @param size number of bootstrap samples.
se.correct <- function(fit, calls, size){
  coefs <- coef(fit)
  npars <- length(coefs)
  res <- matrix(0, nrow = size, ncol = npars + 1)
  traps <- fit$traps
  mask <- fit$mask
  colnames(res) <- c(names(coefs), "maxgrad")
  for (i in 1:size){
    capt <- sim.capt(fit, calls = calls)
    bootfit <- admbsecr(capt, traps = traps, mask = mask)
    res[i, ] <- c(coef(bootfit), bootfit$maxgrad)
  }
  conv <- res[, "maxgrad"] > -1
  res <- res[conv, -which(colnames(res) == "maxgrad")]
  bias <- apply(res, 2, mean) - coefs
  coefficients.corrected <- coefs - bias
  se.corrected <- apply(res, 2, sd)
  out <- list(boots = res, coefficients.corrected = coefficients.corrected,
              se.corrected = se.corrected)
  fit$se.correct <- out
  class(fit) <- c("secorrect", class(fit))
  fit
}

## Modified from R2admb.
#' @S3method summary secorrect
summary.secorrect <- function(object, ...){
  coef.p <- unlist(coef(object,"par"))
  s.err <- object$se.correct$se.corrected
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

#' @S3method stdEr secorrect
stdEr.secorrect <- function(object, ...){
  object$se.correct$se.corrected
}
