#' Bootstrapping SECR data
#'
#' Carries out a parametric bootstrap, based on a model fitted using
#' \link[admbsecr]{admbsecr}. For fits based on acoustic surveys where
#' the argument \code{call.freqs} is provided to the \code{admbsecr}
#' function, the simulated data allocates multiple calls to the same
#' location based on an estimated distribution of the call
#' frequencies. Using a parametric bootstrap is currently the only way
#' parameter uncertainty can be estimated for such models (see
#' Details).
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
    res
}
