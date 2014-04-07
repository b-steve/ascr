#' Extract admbsecr model coefficients
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{admbsecr}.
#'
#' @param object A fitted model from \link[admbsecr]{admbsecr}.
#' @param pars A character vector containing either parameter names,
#' or a subset of \code{"all"}, \code{"derived"}, \code{"fitted"}, and
#' \code{"linked"}; \code{"fitted"} corresponds to the parameters of
#' interest, \code{"derived"} corresponds to quantities that are
#' functions of these parameters (e.g., the effective survey area or
#' animal density from an acoustic survey), and \code{"linked"}
#' corresponds to the parameters AD Model Builder has maximised the
#' likelihood over.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @examples
#' coef(simple.hn.fit)
#' coef(simple.hn.fit, pars = "all")
#' coef(simple.hn.fit, pars = "derived")
#'
#' @method coef admbsecr
#' @S3method coef admbsecr
#'
#' @export
coef.admbsecr <- function(object, pars = "fitted", ...){
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    par.names <- names(object$coefficients)
    if (!all(pars %in% c("fitted", "derived", "linked", par.names))){
        stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
    }
    if (any(c("fitted", "derived", "linked") %in% pars)){
        which.linked <- grep("_link", par.names)
        linked <- object$coefficients[which.linked]
        which.derived <- which(par.names == "esa" | par.names == "Da")
        derived <- object$coefficients[which.derived]
        fitted <- object$coefficients[-c(which.linked, which.derived)]
        out <- mget(pars)
        names(out) <- NULL
        out <- c(out, recursive = TRUE)
    } else {
        out <- object$coefficients[pars]
    }
    out
}

#' @rdname coef.admbsecr
#'
#' @param correct.bias Logical, if \code{TRUE}, estimated biases are
#' subtracted from estimated parameter values.
#'
#' @method coef admbsecr.boot
#' @S3method coef admbsecr.boot
#'
#' @export
coef.admbsecr.boot <- function(object, pars = "fitted",
                               correct.bias = FALSE, ...){
    out <- coef.admbsecr(object, pars)
    if (correct.bias){
        out <- out - get.bias(object, pars)
    }
    out
}
#' Extract the variance-covariance matrix from an admbsecr model
#' object
#'
#' Extracts the variance-covariance matrix for parameters in a model
#' fitted using \link[admbsecr]{admbsecr}.
#'
#' @inheritParams coef.admbsecr
#'
#' @examples
#' vcov(simple.hn.fit)
#' vcov(simple.hn.fit, pars = "all")
#' vcov(simple.hn.fit, pars = "derived")
#'
#' @method vcov admbsecr
#' @S3method vcov admbsecr
#'
#' @export
vcov.admbsecr <- function(object, pars = "fitted", ...){
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    par.names <- names(object$coefficients)
    if (!all(pars %in% c("fitted", "derived", "linked", par.names))){
        stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
    }
    if (any(c("fitted", "derived", "linked") %in% pars)){
        which.linked <- grep("_link", par.names)
        which.derived <- which(par.names == "esa" | par.names == "Da")
        which.fitted <- (1:length(par.names))[-c(which.linked, which.derived)]
        keep <- NULL
        if ("fitted" %in% pars){
            keep <- c(keep, which.fitted)
        }
        if ("derived" %in% pars){
            keep <- c(keep, which.derived)
        }
        if ("linked" %in% pars){
            keep <- c(keep, which.linked)
        }
    } else {
        keep <- pars
    }
    object$vcov[keep, keep, drop = FALSE]
}

#' Extract the variance-covariance matrix from a bootstrapped admbsecr
#' model object
#'
#' Extracts the variance-covariance matrix for parameters in a model
#' fitted using \link[admbsecr]{admbsecr}, with a bootstrap procedure
#' carried out using \link[admbsecr]{boot.admbsecr}.
#'
#' @inheritParams coef.admbsecr
#'
#' @method vcov admbsecr.boot
#' @S3method vcov admbsecr.boot
#'
#' @export
vcov.admbsecr.boot <- function(object, pars = "fitted", ...){
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    par.names <- names(object$coefficients)
    if (!all(pars %in% c("fitted", "derived", "linked", par.names))){
        stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
    }
    if (any(c("fitted", "derived", "linked") %in% pars)){
        which.linked <- grep("_link", par.names)
        which.derived <- which(par.names == "esa" | par.names == "Da")
        which.fitted <- (1:length(par.names))[-c(which.linked, which.derived)]
        keep <- NULL
        if ("fitted" %in% pars){
            keep <- c(keep, which.fitted)
        }
        if ("derived" %in% pars){
            keep <- c(keep, which.derived)
        }
        if ("linked" %in% pars){
            keep <- c(keep, which.linked)
        }
    } else {
        keep <- pars
    }
    object$boot$vcov[keep, keep, drop = FALSE]
}

#' Extract standard errors from an admbsecr model fit
#'
#' Extracts standard errors for estimated and derived parameters from
#' a model fitted using \link[admbsecr]{admbsecr}.
#'
#' @inheritParams coef.admbsecr
#'
#' @examples
#' stdEr(simple.hn.fit)
#' stdEr(simple.hn.fit, pars = "all")
#' stdEr(simple.hn.fit, pars = "derived")
#'
#' @method stdEr admbsecr
#' @S3method stdEr admbsecr
#'
#' @export
stdEr.admbsecr <- function(object, pars = "fitted", ...){
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    par.names <- names(object$coefficients)
    if (!all(pars %in% c("fitted", "derived", "linked", par.names))){
        stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
    }
    if (any(c("fitted", "derived", "linked") %in% pars)){
        which.linked <- grep("_link", par.names)
        linked <- object$se[which.linked]
        which.derived <- which(par.names == "esa" | par.names == "Da")
        derived <- object$se[which.derived]
        fitted <- object$se[-c(which.linked, which.derived)]
        out <- mget(pars)
        names(out) <- NULL
        out <- c(out, recursive = TRUE)
    } else {
        out <- object$se[pars]
    }
    out
}

#' Extract standard errors from a bootstrapped admbsecr model object
#'
#' Extracts standard errors for parameters in a model fitted using
#' \link[admbsecr]{admbsecr}, with a bootstrap procedure carried out
#' using \link[admbsecr]{boot.admbsecr}.
#'
#' @inheritParams coef.admbsecr
#'
#' @method stdEr admbsecr.boot
#' @S3method stdEr admbsecr.boot
#'
#' @export
stdEr.admbsecr.boot <- function(object, pars = "fitted", ...){
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    par.names <- names(object$coefficients)
    if (!all(pars %in% c("fitted", "derived", "linked", par.names))){
        stop("Argument 'pars' must either contain a vector of parameter names, or a subset of \"fitted\", \"derived\", \"linked\", and \"all\".")
    }
    if (any(c("fitted", "derived", "linked") %in% pars)){
        which.linked <- grep("_link", par.names)
        linked <- object$boot$se[which.linked]
        which.derived <- which(par.names == "esa" | par.names == "Da")
        derived <- object$boot$se[which.derived]
        fitted <- object$boot$se[-c(which.linked, which.derived)]
        out <- mget(pars)
        names(out) <- NULL
        out <- c(out, recursive = TRUE)
    } else {
        out <- object$boot$se[pars]
    }
    out
}

#' Extract AIC from an admbsecr model object
#'
#' Extracts the AIC from an admbsecr model object.
#'
#' If the model is based on an acoustic survey where there are
#' multiple calls per individual, then AIC should not be used for
#' model selection. This function therefore returns NA in this case.
#'
#' @inheritParams coef.admbsecr
#' @inheritParams stats::AIC
#'
#' @method AIC admbsecr
#' @S3method AIC admbsecr
#'
#' @export
AIC.admbsecr <- function(object, ..., k = 2){
    if (object$fit.freqs){
        out <- NA
    } else {
        out <- deviance(object) + k*length(coef(object))
    }
    out
}

#' Summarising admbsecr model fits
#'
#' Provides a useful summary of the model fit.
#'
#' @inheritParams coef.admbsecr
#'
#' @method summary admbsecr
#' @S3method summary admbsecr
#'
#' @export
summary.admbsecr <- function(object, ...){
    coefs <- coef(object, "fitted")
    derived <- coef(object, "derived")
    coefs.se <- stdEr(object, "fitted")
    derived.se <- stdEr(object, "derived")
    out <- list(coefs = coefs, derived = derived, coefs.se = coefs.se,
                derived.se = derived.se)
    class(out) <- c("summary.admbsecr", class(out))
    out
}

#' @method print summary.admbsecr
#' @S3method print summary.admbsecr
print.summary.admbsecr <- function(x, ...){
    n.coefs <- length(x$coefs)
    n.derived <- length(x$derived)
    mat <- matrix(0, nrow = n.coefs + n.derived + 1, ncol = 2)
    mat[1:n.coefs, 1] <- c(x$coefs)
    mat[1:n.coefs, 2] <- c(x$coefs.se)
    mat[n.coefs + 1, ] <- NA
    mat[(n.coefs + 2):(n.coefs + n.derived + 1), ] <- c(x$derived, x$derived.se)
    rownames(mat) <- c(names(x$coefs), "---", names(x$derived))
    colnames(mat) <- c("Estimate", "Std. Error")
    cat("Coefficients:", "\n")
    printCoefmat(mat, na.print = "")
}

#' Confidence intervals for admbsecr model parameters
#'
#' Computes confidence intervals for one or more parameters estimated
#' in an admbsecr model object.
#'
#' Options for the argument \code{method} are as follows:
#' \code{"default"} for intervals based on a normal approximation
#' using the calculated standard errors (for objects of class
#' \code{admbsecr.boot}, these standard errors are calculated from the
#' bootstrap procedure); \code{"default.bc"} is a bias-corrected
#' version of \code{default}, whereby the estimated bias is subtracted
#' from each confidence limit; \code{"basic"} for the so-called
#' "basic" bootstrap method; and \code{"percentile"} for intervals
#' calculated using the bootstrap percentile method (the latter three
#' are only available for objects of class \code{admbsecr.boot}; see
#' Davison & Hinkley, 1997, for details).
#'
#' For method \code{"default"} with objects of class
#' \code{admbsecr.boot}, the appropriateness of the normal
#' approximation can be evaluated by setting \code{qqnorm} to
#' \code{TRUE}. If this indicates a poor fit, set \code{linked} to
#' \code{TRUE} and evaluate the QQ plot to see if this yields an
#' improvement (see Davison & Hinkley, 1997, pp. 194, for details).
#'
#' @references Davison, A. C., and Hinkley, D. V. (1997)
#' \emph{Bootstrap methods and their application}. Cambridge:
#' Cambridge University Press.
#'
#' @param parm A character vector containing either parameter names,
#' or a subset of \code{"all"}, \code{"derived"}, \code{"fitted"}, and
#' \code{"linked"}; \code{"fitted"} corresponds to the parameters of
#' interest, \code{"derived"} corresponds to quantities that are
#' functions of these parameters (e.g., the effective survey area or
#' animal density from an acoustic survey), and \code{"linked"}
#' corresponds to the parameters AD Model Builder has maximised the
#' likelihood over.
#' @param linked Logical, if \code{TRUE}, intervals for fitted
#' parameters are calculated on their link scales, then transformed
#' back onto their "real" scales.
#' @inheritParams coef.admbsecr
#' @inheritParams stats::confint
#'
#' @method confint admbsecr
#' @S3method confint admbsecr
#'
#' @export
confint.admbsecr <- function(object, parm = "fitted", level = 0.95, linked = FALSE, ...){
    if (!object$args$hess){
        stop("Standard errors not calculated; use boot.admbsecr() or refit with 'hess = TRUE', if appropriate.")
    }
    calc.cis(object, parm, level, method = "default", linked, qqplot = FALSE,
             boot = FALSE, ask = FALSE, ...)
}

#' @param method A character string specifying the method used to
#' calculate the confidence intervals. See 'Details' below.
#' @param qqplot Logical, if \code{TRUE} and \code{method} is
#' \code{"default"} then a normal QQ plot is plotted. The default
#' method is based on a normal approximation; this plot tests its
#' validity.
#' @param ask Logical, if \code{TRUE}, hitting return will show the
#' next plot.
#'
#' @rdname confint.admbsecr
#' @method confint admbsecr.boot
#' @S3method confint admbsecr.boot
#'
#' @export
confint.admbsecr.boot <- function(object, parm = "fitted", level = 0.95, method = "default",
                                  linked = FALSE, qqplot = FALSE, ask = TRUE, ...){
    calc.cis(object, parm, level, method, linked, qqplot, boot = TRUE, ask, ...)
}

calc.cis <- function(object, parm, level, method, linked, qqplot, boot, ask, ...){
    if (any(c("all", "derived", "fitted") %in% parm)){
        parm <- names(coef(object, pars = parm))
    }
    if (linked){
        fitted.names <- names(coef(object, "fitted")[parm])
        fitted.names <- fitted.names[fitted.names != "mu.freqs"]
        linked.names <- paste(fitted.names, "_link", sep = "")
        link.parm <- linked.names[!(linked.names %in% parm)]
        all.parm <- c(parm, link.parm)
    } else {
        all.parm <- parm
    }
    if (method == "default" | method == "default.bc"){
        mat <- cbind(coef(object, pars = "all")[all.parm],
                     stdEr(object, pars = "all")[all.parm])
        FUN.default <- function(x, level){
            x[1] + qnorm((1 - level)/2)*c(1, -1)*x[2]
        }
        out <- t(apply(mat, 1, FUN.default, level = level))
        if (method == "default.bc"){
            out <- out - get.bias(object, parm)
        }
        if (qqplot & boot){
            opar <- par(ask = ask)
            for (i in parm){
                if (linked){
                    if (i %in% fitted.names){
                        j <- linked.names[fitted.names == i]
                    }
                } else {
                    j <- i
                }
                qqnorm(object$boot$boots[, j], main = i)
                abline(mean(object$boot$boots[, j], na.rm = TRUE),
                       sd(object$boot$boots[, j], na.rm = TRUE))
            }
            par(opar)
        }
    } else if (method == "basic"){
        qs <- t(apply(object$boot$boots[, all.parm, drop = FALSE], 2, quantile,
                      probs = c((1 - level)/2, 1 - (1 - level)/2),
                      na.rm = TRUE))
        mat <- cbind(coef(object, pars = "all")[all.parm], qs)
        FUN.basic <- function(x){
            2*x[1] - c(x[3], x[2])
        }
        out <- t(apply(mat, 1, FUN.basic))
    } else if (method == "percentile"){
        out <- t(apply(object$boot$boots[, all.parm, drop = FALSE], 2, quantile,
                       probs = c((1 - level)/2, 1 - (1 - level)/2),
                       na.rm = TRUE))
    }
    if (linked){
        for (i in fitted.names){
            linked.name <- paste(i, "_link", sep = "")
            out[i, ] <- object$par.unlinks[[i]](out[linked.name, ])
        }
        out <- out[parm, , drop = FALSE]
    }
    percs <- c(100*(1 - level)/2, 100*(1 - (1 - level)/2))
    colnames(out) <- paste(round(percs, 2), "%")
    out
}
