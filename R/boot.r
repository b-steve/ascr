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
#' Note that generic functions \link{stdEr} and \link{vcov} with an
#' object of class \code{admbsecr.boot} as the main argument will
#' return standard errors and the variance-covariance matrix for
#' estimated parameters \emph{based on the bootstrap procedure} (via
#' the \link{stdEr.admbsecr.boot} and \link{vcov.admbsecr.boot}
#' methods). For standard errors and the variance-covariance matrix
#' based on maximum likelihood asymptotic theory, the methods
#' \link{stdEr.admbsecr} and \link{vcov.admbsecr} must be called
#' directly.
#'
#' @section Bootstrapping for acoustic surveys:
#'
#' For fits based on acoustic surveys where the argument
#' \code{call.freqs} is provided to the \code{admbsecr} function, the
#' simulated data allocates multiple calls to the same location based
#' on an estimated distribution of the call frequencies. Using a
#' parametric bootstrap is currently the only way parameter
#' uncertainty can be estimated for such models (see Stevenson et al.,
#' in prep., for details).
#'
#' @section Monte Carlo error:
#'
#' There will be some error in esimates based on the parametric
#' bootstrap (e.g., standard errors and estimates of bias) because the
#' number of bootstrap simulations is not infinite. By default, this
#' function calculates Monte Carlo error using a bootstrap of the
#' estimates obtained from the initial bootstrap procedure; see
#' Equation (9) in Koehler, Brown and Haneuse (2009). Note that this
#' secondary bootstrap does not require the fitting of any further
#' models, and so the increased processing time due to this procedure
#' is negligible.
#'
#' Monte Carlo error for standard errors and estimates of bias can be
#' extracted using the function \link{get.mce}.
#'
#' @references Koehler, E., Brown, E., and Haneuse, S. J.-P. A. (2009)
#' On the assessment of Monte Carlo error in sumulation-based
#' statistical analyses. \emph{The American Statistician},
#' \strong{63}: 155--162.
#'
#' @references Stevenson, B. C., Borchers, D. L., Altwegg, R., Measey,
#' G. J., Swift, R. J., and Gillespie, D. M. (in prep.) An acoustic
#' spatially explicit capture-recapture method for estimating
#' vocalizing amphibian density.
#'
#' @return A list of class \code{"admbsecr.boot"}. Components contain
#' information such as estimated parameters and standard errors. The
#' best way to access such information, however, is through the
#' variety of helper functions provided by the admbsecr package. S3
#' methods \link{stdEr.admbsecr.boot} and \link{vcov.admbsecr.boot}
#' can be used to return standard errors and the variance-covariance
#' matrix of estimated parameters based on the bootstrap procedure.
#'
#' @param fit A fitted \code{admbsecr} model object.
#' @param N The number of bootstrap resamples.
#' @param prog Logical, if \code{TRUE}, a progress bar is shown. Only
#' available if \code{n.cores} is 1.
#' @param n.cores A positive integer representing the number of cores
#' to use for parallel processing.
#' @param M The number of bootstrap resamples for the secondary
#' bootstrap used to calculate Monte Carlo error. See 'Details'
#' below. If M = 0, then this is skipped.
#'
#'
#' @examples
#' \dontrun{
#' ## In practice, N should be >> 100, but this leads to long computation time for a simple example.
#' boot.fit <- boot.admbsecr(fit = simple.hn.fit, N = 100)
#' }
#'
#' @export
boot.admbsecr <- function(fit, N, prog = TRUE, n.cores = 1, M = 10000){
    args <- fit$args
    orig.sv <- args$sv
    ## Set start values to estimated parameters.
    args$sv <- get.par(fit, "fitted", as.list = TRUE)
    ## Removing scalefactors.
    args$sf <- NULL
    ## Setting trace to false.
    args$trace <- FALSE
    ## Don't calculate the Hessian for bootstrap fits.
    args$hess <- FALSE
    ## Setting start value for g0 away from 1.
    if ("g0" %in% names(args$sv)){
        args$sv[["g0"]] <- min(c(0.95, args$sv[["g0"]]))
    }
    call.freqs <- args$call.freqs
    coefs <- fit$coefficients
    par.names <- names(coefs)
    n.pars <- length(coefs)
    seeds <- sample(1:1e8, size = N)
    seed.mce <- sample(1:1e8, size = 1)
    ## Function to get fit.boot.
    FUN <- function(i, fit, args, call.freqs, seeds, prog){
        set.seed(seeds[i])
        ## Simulating capture history.
        args$capt <- sim.capt(fit)
        ## Simulating calling frequencies (if required).
        if (fit$fit.freqs){
            if (length(call.freqs) > 1){
                args$call.freqs <- sample(call.freqs, replace = TRUE)
            }
        }
        ## Fitting model.
        fit.boot <- try(do.call("admbsecr", args), silent = TRUE)
        ## If unconverged, refit model with default start values.
        if (fit.boot$maxgrad < -0.01 | "try-error" %in% class(fit.boot)){
            args$sv <- NULL
            fit.boot <- try(do.call("admbsecr", args), silent = TRUE)
        }
        ## If still unconverged, give up and report NA.
        if (fit.boot$maxgrad < -0.01 | "try-error" %in% class(fit.boot)){
            n.par <- length(fit$coefficients)
            out <- rep(NA, n.par + 1)
        } else {
            out <- c(fit.boot$coefficients, fit.boot$maxgrad)
        }
        if (prog){
            cat(i, "\n", file = "prog.txt", append = TRUE)
        }
        out
    }
    if (n.cores == 1){
        res <- matrix(0, nrow = N, ncol = n.pars + 1)
        colnames(res) <- c(par.names, "maxgrad")
        ## Setting up progress bar.
        if (prog){
            pb <- txtProgressBar(min = 0, max = N, style = 3)
        }
        for (i in 1:N){
            res[i, ] <- FUN(i, fit = fit, args = args, call.freqs = call.freqs,
                            seeds = seeds, prog = FALSE)
            ## Updating progress bar.
            if (prog){
                setTxtProgressBar(pb, i)
            }
        }
        ## Closing progress bar.
        if (prog){
            close(pb)
        }
    } else {
        if (!require(parallel)){
            stop("The parallel package is required for n.cores > 1. Please install.")
        }
        if (n.cores > detectCores()){
            stop("The argument n.cores is greater than the number of available cores.")
        }
        cluster <- makeCluster(n.cores)
        clusterEvalQ(cluster, {
            library(admbsecr)
        })
        if (prog){
            file.create("prog.txt")
        }
        res <- t(parSapplyLB(cluster, 1:N, FUN, fit = fit, args = args,
                             call.freqs = call.freqs, seeds = seeds, prog = prog))
        stopCluster(cluster)
        if (prog){
            unlink("prog.txt")
        }
    }
    ## Calculating bootstrapped standard errors, correlations and
    ## covariances.
    maxgrads <- res[, ncol(res)]
    res <- res[, -ncol(res)]
    se <- apply(res, 2, sd, na.rm = TRUE)
    names(se) <- par.names
    cor <- diag(n.pars)
    dimnames(cor) <- list(par.names, par.names)
    vcov <- diag(se^2)
    dimnames(vcov) <- list(par.names, par.names)
    for (i in 1:(n.pars - 1)){
        for (j in (i + 1):n.pars){
            cor[i, j] <- cor[j, i] <- cor(res[, i], res[, j], use = "na.or.complete")
            vcov[i, j] <- vcov[j, i] <- cor[i, j]*se[i]*se[j]
        }
    }
    bias <- apply(res, 2, mean, na.rm = TRUE) - coefs
    ## Bootstrap to calculate MCE for bias and standard errors.
    bias.mce <- se.mce <- numeric(n.pars)
    names(bias.mce) <- names(se.mce) <- par.names
    if (M > 0){
        set.seed(seed.mce)
        converged <- which(!is.na(res[, 1]))
        n.converged <- length(converged)
        mce.boot <- matrix(sample(converged, size = n.converged*M,
                                  replace = TRUE), nrow = M,
                           ncol = n.converged)
        for (i in par.names){
            par.boot <- matrix(res[mce.boot, i], nrow = M, ncol = n.converged)
            bias.mce[i] <- sd(apply(par.boot, 1, mean) - coefs[i] - bias[i])
            se.mce[i] <- sd(apply(par.boot, 1, sd))
        }
    } else {
        bias.mce <- NA
        se.mce <- NA
    }
    out <- fit
    boot <- list(boots = res, se = se, se.mce = se.mce, cor = cor, vcov = vcov,
                 bias = bias, bias.mce = bias.mce, maxgrads = maxgrads)
    out$boot <- boot
    class(out) <- c("admbsecr.boot", class(fit))
    out
}
