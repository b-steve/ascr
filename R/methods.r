#' Extract admbsecr model coefficients
#'
#' Extracts estimated and derived coefficients from a model fitted
#' using \link[admbsecr]{admbsecr}.
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
#' @method coef admbsecr
#' @S3method coef admbsecr
coef.admbsecr <- function(object, pars = "fitted", ...){
    if (pars == "all"){
        out <- object$coefficients
    } else if (pars == "fitted"){
        out <- object$coefficients[object$coefficients != "esa"]
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
#' @method vcov admbsecr
#' @S3method vcov admbsecr
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
