#' Bootstrapping SECR data
#'
#' Carries out a parametric bootstrap, based on a model fitted using
#' \link{admbsecr}.
#'
#' For each bootstrap resample, a new population of individuals is
#' simulated within the mask area. Detections of these individuals are
#' simulated using the estimated detection function. For detected
#' individuals, additional informatin is simulated from the estimated
#' distribution of measurement error. The original model is then
#' re-fitted to these simulated data, and parameter estimates for each
#' iteration saved in the component \code{boot} of the returned list.
#' 
#' For fits based on acoustic surveys where the argument
#' \code{call.freqs} is provided to the \code{admbsecr} function, the
#' simulated data allocates multiple calls to the same location based
#' on an estimated distribution of the call frequencies. Using a
#' parametric bootstrap is currently the only way parameter
#' uncertainty can be estimated for such models.
#'
#' @return A list of class \code{"admbsecr.boot"}. Components contain
#' information such as estimated parameters and standard errors. The
#' best way to access such information, however, is through the
#' variety of helper functions provided by the admbsecr package. S3
#' methods \link{stdEr.admbsecr.boot} and \link{vcov.admbsecr.boot}
#' can be used to return values of based on the bootstrap procedure.
#' 
#' @param fit A fitted \code{admbsecr} model object.
#' @param N The number of bootstrap resamples.
#' @param prog Logical, if \code{TRUE}, a progress bar is shown.
#' @export
boot.admbsecr <- function(fit, N, prog = TRUE){
    args <- fit$args
    orig.sv <- args$sv
    ## Set start values to estimated parameters.
    args$sv <- get.par(fit, "all", as.list = TRUE)
    ## Logical value for call frequency data.
    fit.freqs <- fit$fit.freqs
    if (fit.freqs){
        call.freqs <- args$call.freqs
    }
    n.pars <- length(fit$coefficients)
    res <- matrix(0, nrow = N, ncol = n.pars)
    colnames(res) <- names(fit$coefficients)
    ## Setting up progress bar.
    if (prog){
        pb <- txtProgressBar(min = 0, max = N, style = 3)
    }
    for (i in 1:N){
        ## Simulating capture history.
        args$capt <- sim.capt(fit)
        ## Simulating calling frequencies (if required).
        if (fit.freqs){
            if (length(call.freqs) > 1){
                args$call.freqs <- sample(call.freqs, replace = TRUE)
            }
        }
        ## Fitting model.
        fit.boot <- do.call("admbsecr", args)
        if (fit.boot$maxgrad < -0.01){
            res[i, ] <- NA
        } else {
            res[i, ] <- fit.boot$coefficients
        }
        ## Updating progress bar.
        if (prog){
            setTxtProgressBar(pb, i)
        }
    }
    ## Closing progress bar.
    if (prog){
        close(pb)
    }
    ## Calculating bootstrapped standard errors, correlations and
    ## covariances.
    se <- apply(res, 2, sd)
    names(se) <- names(fit$se)
    cor <- diag(n.pars)
    dimnames(cor) <- dimnames(fit$cor)
    vcov <- diag(se^2)
    dimnames(vcov) <- dimnames(fit$vcov)
    for (i in 1:(n.pars - 1)){
        for (j in (i + 1):n.pars){
            cor[i, j] <- cor[j, i] <- cor(res[, i], res[, j])
            vcov[i, j] <- vcov[j, i] <- cor[i, j]*se[i]*se[j]
        }
    }
    out <- fit
    out$boot.se <- se
    out$boot.cor <- cor
    out$boot.vcov <- vcov
    out$boot <- res
    class(out) <- c("admbsecr.boot", class(fit))
    out
}
