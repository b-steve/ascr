#' Extract admbsecr model coefficients
#'
#' Extracts estimated and derived parameters from a model fitted using
#' \link{admbsecr}.
#'
#' @param object A fitted model from \link[admbsecr]{admbsecr}.
#' @param pars A character containing a subset of \code{"all"},
#' \code{"derived"}, \code{"fitted"}, and \code{"linked"};
#' \code{"fitted"} corresponds to the parameters of interest,
#' \code{"derived"} corresponds to quantities that are functions of
#' these parameters (e.g., the effective survey area), and
#' \code{"linked"} corresponds to the parameters AD Model Builder has
#' maximised the likelihood over.
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
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    if (!all(pars %in% c("fitted", "derived", "linked"))){
        stop("Argument 'pars' must contain a subset of \"fitted\", \"derived\", and \"linked\"")
    }
    par.names <- names(object$coefficients)
    which.linked <- grep("_link", par.names)
    linked <- object$coefficients[which.linked]
    which.derived <- which(par.names == "esa")
    derived <- object$coefficients[which.derived]
    fitted <- object$coefficients[-c(which.linked, which.derived)]
    out <- mget(pars)
    names(out) <- NULL
    c(out, recursive = TRUE)
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
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    if (!all(pars %in% c("fitted", "derived", "linked"))){
        stop("Argument 'pars' must contain a subset of \"fitted\", \"derived\", and \"linked\"")
    }
    par.names <- names(object$coefficients)
    keep <- NULL
    which.linked <- grep("_link", par.names)
    which.derived <- which(par.names == "esa")
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
    object$vcov[keep, keep]
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
    if ("all" %in% pars){
        pars <- c("fitted", "derived", "linked")
    }
    if (!all(pars %in% c("fitted", "derived", "linked"))){
        stop("Argument 'pars' must contain a subset of \"fitted\", \"derived\", and \"linked\"")
    }
    par.names <- names(object$coefficients)
    which.linked <- grep("_link", par.names)
    linked <- object$se[which.linked]
    which.derived <- which(par.names == "esa")
    derived <- object$se[which.derived]
    fitted <- object$se[-c(which.linked, which.derived)]
    out <- mget(pars)
    names(out) <- NULL
    c(out, recursive = TRUE)
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
