#' Extract admbsecr model coefficients
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{admbsecr}.
#'
#' @param object A fitted model from \link[admbsecr]{admbsecr}.
#' @param pars A character string; either \code{"all"},
#' \code{"fitted"}, or \code{"derived"}. If \code{"all"}, all
#' estimated parameters are returned. If \code{"fitted"}, only those
#' parameters for which the likelihood was maximised over are
#' returned. If \code{"derived"}, only those parameters (e.g., the
#' effective survey area) which are functions of fitted parameters are
#' returned.
#' @param ... Other parameters (for S3 generic compatibility).
#'
#' @examples
#' coef(simple.hn.fit)
#' coef(simple.hn.fit, pars = "all")
#' coef(simple.hn.fit, pars = "derived")
#' 
#' @method coef admbsecr
#' @S3method coef admbsecr
#' @export
coef.admbsecr <- function(object, pars = "fitted", ...){
    if (pars == "all"){
        out <- object$coefficients
    } else if (pars == "fitted"){
        out <- object$coefficients[names(object$coefficients) != "esa"]
    } else if (pars == "derived"){
        out <- object$coefficients["esa"]
    } else {
        stop("Argument 'pars' must be either \"all\", \"fitted\", or \"derived\".")
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
#' @export
vcov.admbsecr <- function(object, pars = "fitted", ...){
    if (pars == "all"){
        out <- object$vcov
    } else if (pars == "fitted"){
        par.names <- names(object$coefficients)
        keep.pars <- par.names[par.names != "esa"]
        out <- object$vcov[keep.pars, keep.pars]
    } else if (pars == "derived"){
        out <- object$vcov["esa", "esa"]
    } else {
        stop("Argument 'pars' must be either \"all\", \"fitted\", or \"derived\".")
    }
    out
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
vcov.admbsecr.boot <- function(object, pars = "fitted", ...){
    if (pars == "all"){
        out <- object$boot.vcov
    } else if (pars == "fitted"){
        par.names <- names(object$coefficients)
        keep.pars <- par.names[par.names != "esa"]
        out <- object$boot.vcov[keep.pars, keep.pars]
    } else if (pars == "derived"){
        out <- object$boot.vcov["esa", "esa"]
    } else {
        stop("Argument 'pars' must be either \"all\", \"fitted\", or \"derived\".")
    }
    out
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
#' @export
stdEr.admbsecr <- function(object, pars = "fitted", ...){
    if (pars == "all"){
        out <- object$se
    } else if (pars == "fitted"){
        out <- object$se[names(object$se) != "esa"]
    } else if (pars == "derived"){
        out <- object$se["esa"]
    } else {
        stop("Argument 'pars' must be either \"all\", \"fitted\", or \"derived\".")
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
stdEr.admbsecr.boot <- function(object, pars = "fitted", ...){
    if (pars == "all"){
        out <- object$boot.se
    } else if (pars == "fitted"){
        out <- object$boot.se[names(object$se) != "esa"]
    } else if (pars == "derived"){
        out <- object$boot.se["esa"]
    } else {
        stop("Argument 'pars' must be either \"all\", \"fitted\", or \"derived\".")
    }
    out
}
