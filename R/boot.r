#' Bootstrapping SECR data
#'
#' Carries out a parametric bootstrap, based on a model fitted using
#' \link{fit.ascr}.
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
#' object of class \code{ascr.boot} as the main argument will
#' return standard errors and the variance-covariance matrix for
#' estimated parameters \emph{based on the bootstrap procedure} (via
#' the \link{stdEr.ascr.boot} and \link{vcov.ascr.boot}
#' methods). For standard errors and the variance-covariance matrix
#' based on maximum likelihood asymptotic theory, the methods
#' \link{stdEr.ascr} and \link{vcov.ascr} must be called
#' directly.
#'
#' If \code{infotypes} is provided it should take the form of a list,
#' where each component is a subset of information types (i.e.,
#' \code{fit$infotypes}) used to fit the original model. A \code{NULL}
#' component is associated with no additional information. The
#' bootstrap procedure is then repeated for each component, only
#' utilising the appropriate additional information. In practice this
#' is only useful if the user is looking to investigate the benefits
#' of including particular information types. The results from these
#' extra bootstrap procedures can be found in the
#' \code{boot$extra.boots} component of the returned object.
#'
#' the bootstrap procedure will be
#' repeated for each subset; usually this is only useful to
#' investigate the impact on the density estimator associated with
#' various combinations of information types. A \code{NULL} component
#' will bootstrap without any additional information types.
#'
#' @section Bootstrapping for acoustic surveys:
#'
#' For fits based on acoustic surveys where the argument
#' \code{cue.rates} is provided to the \code{fit.ascr} function, the
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
#'     On the assessment of Monte Carlo error in sumulation-based
#'     statistical analyses. \emph{The American Statistician},
#'     \strong{63}: 155--162.
#'
#' @references Stevenson, B. C., Borchers, D. L., Altwegg, R., Measey,
#'     G. J., Swift, R. J., and Gillespie, D. M. (in prep.) A general
#'     framework for animal density estimation from acoustic detection
#'     data.
#'
#' @return A list of class \code{"ascr.boot"}. Components contain
#'     information such as estimated parameters and standard
#'     errors. The best way to access such information, however, is
#'     through the variety of helper functions provided by the
#'     ascr package. S3 methods \link{stdEr.ascr.boot} and
#'     \link{vcov.ascr.boot} can be used to return standard errors
#'     and the variance-covariance matrix of estimated parameters
#'     based on the bootstrap procedure.
#'
#' @param fit A fitted \code{ascr} model object.
#' @param N The number of bootstrap resamples.
#' @param prog Logical, if \code{TRUE}, a progress bar is shown. Only
#'     available if \code{n.cores} is 1.
#' @param n.cores A positive integer representing the number of cores
#'     to use for parallel processing.
#' @param M The number of bootstrap resamples for the secondary
#'     bootstrap used to calculate Monte Carlo error. See 'Details'
#'     below. If M = 0, then this is skipped.
#' @param infotypes A list, where each component contains information
#'     types for subsequent bootstrap procedures. See 'Details'.
#'
#'
#'
#' @examples
#' \dontrun{
#' ## In practice, N should be >> 100, but this leads to long computation time for a simple example.
#' boot.fit <- boot.ascr(fit = example$fits$simple.hn, N = 100)
#' }
#'
#' @export
boot.ascr <- function(fit, N, prog = TRUE, n.cores = 1, M = 10000, infotypes = NULL){
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
    cue.rates <- args$cue.rates
    coefs <- fit$coefficients
    par.names <- names(coefs)
    n.pars <- length(coefs)
    seeds <- sample(1:1e8, size = N)
    seed.mce <- sample(1:1e8, size = 1)
    ## Function to get fit.boot.
    FUN <- function(i, fit, args, cue.rates, infotypes, seeds, prog){
        set.seed(seeds[i])
        if (fit$n.sessions > 1){
            ## Simulating capture history.
            args$capt <- lapply(sim.capt(fit), function(x) x[c("bincapt", infotypes)])
            n.dets <- sum(sapply(args$capt, function(x) nrow(x$bincapt)))
        } else {
            ## Simulating capture history.
            args$capt <- sim.capt(fit)[c("bincapt", infotypes)]
            n.dets <- nrow(args$capt$bincapt)
        }
        ## If no calls simulated, set density to 0 and other parameters to NA.
        if (n.dets == 0){
            n.par <- length(fit$coefficients)
            out <- rep(NA, n.par + 1)
            out[names(fit$coefficients) %in% c("D", "Da", "Dc")] <- 0
            out[names(fit$coefficients) == "D_link"] <- -Inf
        } else {
            ## Simulating calling frequencies (if required).
            if (fit$fit.freqs){
                if (length(cue.rates) > 1){
                    args$cue.rates <- sample(cue.rates, replace = TRUE)
                }
            }
            ## Fitting model.
            fit.boot <- suppressWarnings(try(do.call("fit.ascr", args), silent = TRUE))
            ## If unconverged, refit model with default start values.
            if ("try-error" %in% class(fit.boot) || fit.boot$maxgrad < -0.01){
                args$sv <- NULL
                fit.boot <- suppressWarnings(try(do.call("fit.ascr", args), silent = TRUE))
            }
            ## If still unconverged, give up and report NA.
            if ("try-error" %in% class(fit.boot) || fit.boot$maxgrad < -0.01){
                n.par <- length(fit$coefficients)
                out <- rep(NA, n.par + 1)
            } else {
                out <- c(fit.boot$coefficients, fit.boot$maxgrad)
            }
        }
        if (prog){
            cat(i, "\n", file = "prog.txt", append = TRUE)
        }
        out
    }
    if (n.cores == 1){
        ## Main bootstrap.
        res <- matrix(0, nrow = N, ncol = n.pars + 1)
        colnames(res) <- c(par.names, "maxgrad")
        ## Setting up progress bar.
        if (prog){
            pb <- txtProgressBar(min = 0, max = N, style = 3)
        }
        for (i in 1:N){
            res[i, ] <- FUN(i, fit = fit, args = args, cue.rates = cue.rates,
                            infotypes = fit$infotypes, seeds = seeds, prog = FALSE)
            ## Updating progress bar.
            if (prog){
                setTxtProgressBar(pb, i)
            }
        }
        ## Closing progress bar.
        if (prog){
            close(pb)
        }
        ## Additional bootstraps.
        extra.res <- vector(mode = "list", length = length(infotypes))
        names(extra.res) <- names(infotypes)
        for (i in seq(from = 1, by = 1, along.with = infotypes)){
            new.args <- args
            new.args$capt <- args$capt[c("bincapt", infotypes[[i]])]
            new.fit <- suppressWarnings(do.call("fit.ascr", new.args))
            new.n.pars <- length(new.fit$coefficients)
            new.par.names <- names(new.fit$coefficients)
            extra.res[[i]] <- matrix(0, nrow = N, ncol = new.n.pars + 1)
            colnames(extra.res[[i]]) <- c(new.par.names, "maxgrad")
            ## Setting up another progress bar.
            if (prog){
                pb <- txtProgressBar(min = 0, max = N, style = 3)
            }
            for (j in 1:N){
                extra.res[[i]][j, ] <- suppressWarnings(FUN(j, fit = fit, args = args,
                                                            cue.rates = cue.rates,
                                                            infotypes = infotypes[[i]],
                                                            seeds = seeds, prog = FALSE))
                ## Updating progress bar.
                if (prog){
                    setTxtProgressBar(pb, j)
                }
            }
            ## Closing progress bar.
            if (prog){
                close(pb)
            }
        }
    } else {
        if (n.cores > detectCores()){
            stop("The argument n.cores is greater than the number of available cores.")
        }
        cluster <- makeCluster(n.cores)
        clusterEvalQ(cluster, {
            library(ascr)
        })
        ## Main bootstrap.
        if (prog){
            file.create("prog.txt")
        }
        ## IF THERE'S AN ERROR HERE YOU NEED TO REBUILD THE PACKAGE
        ## DUE TO library() CALL ABOVE.
        res <- t(parSapplyLB(cluster, 1:N, FUN, fit = fit, args = args,
                             cue.rates = cue.rates, infotypes = fit$infotypes,
                             seeds = seeds, prog = prog))
        if (prog){
            unlink("prog.txt")
        }
        ## Additional bootstrap.
        extra.res <- vector(mode = "list", length = length(infotypes))
        names(extra.res) <- names(infotypes)
        for (i in seq(from = 1, by = 1, along.with = infotypes)){
            if (prog){
                file.create("prog.txt")
            }
            extra.res[[i]] <- t(parSapplyLB(cluster, 1:N, FUN, fit = fit,
                                            args = args, cue.rates = cue.rates,
                                            infotypes = infotypes[[i]],
                                            seeds = seeds, prog = prog))
            ## Removing maximum gradient component.
            extra.res[[i]] <- extra.res[[i]][, -ncol(extra.res[[i]])]
            if (prog){
                unlink("prog.txt")
            }
        }
        stopCluster(cluster)
    }
    ## Calculating bootstrapped standard errors, correlations and
    ## covariances.
    maxgrads <- res[, ncol(res)]
    ## Removing maximum gradient component.
    res <- res[, -ncol(res), drop = FALSE]
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
                 bias = bias, bias.mce = bias.mce, maxgrads = maxgrads,
                 extra.boots = extra.res)
    out$boot <- boot
    class(out) <- c("ascr.boot", class(fit))
    out
}

## Aliasing old boot.admbsecr() function name.
#' @rdname boot.ascr
#' @export
boot.admbsecr <- boot.ascr

#' Combining subsamples to obtain a standard error.
#'
#' Calculates a single standard error for a parameter that has been
#' calculated by averaging over subsamples.
#'
#' @param ... A number of bootstrap model objects.
#' @param par A character string providing the parameter for which to
#'     calculate a standard error.
#' @param plot Logical, if \code{TRUE}, a boxplot is produced.
#' @param ceiling A threshold value; bootstrapped parameter values
#'     above this are discarded.
subsample.se <- function(..., par, plot = TRUE, ceiling = NULL){
    boot.list <- list(...)
    n.fits <- length(boot.list)
    FUN <- function(x, par){
        x$boot$boots[, par]
    }
    mean.pars <- apply(t(laply(boot.list, FUN, par = par)), 1, mean)
    if (!is.null(ceiling)){
        mean.pars <- mean.pars[mean.pars <= ceiling]
    }
    if (plot){
        boxplot(mean.pars)
    }
    sd(mean.pars, na.rm = TRUE)
}
