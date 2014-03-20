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
#' For fits based on acoustic surveys where the argument
#' \code{call.freqs} is provided to the \code{admbsecr} function, the
#' simulated data allocates multiple calls to the same location based
#' on an estimated distribution of the call frequencies. Using a
#' parametric bootstrap is currently the only way parameter
#' uncertainty can be estimated for such models.
#'
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
#' @param prog Logical, if \code{TRUE}, a progress bar is shown when
#' \code{n.cores} is 1. If \code{n.cores} > 1, then progress is
#' reported in a text file \code{prog.txt} which is created in the
#' working directory.
#' @param n.cores A positive integer representing the number of cores
#' to use for parallel processing.
#'
#' @examples
#' \dontrun{
#' ## In practice, N should be >> 100, but this leads to long computation time for a simple example.
#' boot.fit <- boot.admbsecr(fit = simple.hn.fit, N = 100)
#' }
#' 
#' @export
boot.admbsecr <- function(fit, N, prog = TRUE, n.cores = 1){
    args <- fit$args
    orig.sv <- args$sv
    ## Set start values to estimated parameters.
    args$sv <- get.par(fit, "all", as.list = TRUE)
    call.freqs <- args$call.freqs
    n.pars <- length(fit$coefficients)
    seeds <- sample(1:1e8, size = N)
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
        fit.boot <- do.call("admbsecr", args)
        if (fit.boot$maxgrad < -0.01){
            out <- NA
        } else {
            out <- fit.boot$coefficients
        }
        ## Writing progress bar to progress text file.
        if (prog & file.exists("prog.txt")){
            curr.prog <- strsplit(readLines("prog.txt"), " /")[[1]][1]
            curr.prog <- as.numeric(strsplit(curr.prog, ", ")[[1]][2])
            new.prog <- curr.prog + 1
            N <- length(seeds)
            n.increments <- round(70*new.prog/N)
            perc <- 100*new.prog/N
            cat("  |", rep("=", n.increments), rep(" ", 70 - n.increments),
                "|", rep(" ", 4 - nchar(perc)), perc, "%, ", new.prog, " / ",
                N, "\n", sep = "", file = "prog.txt")
        }
        out
    }  
    if (n.cores == 1){
        res <- matrix(0, nrow = N, ncol = n.pars)
        colnames(res) <- names(fit$coefficients)
        ## Setting up progress bar.
        if (prog){
            pb <- txtProgressBar(min = 0, max = N, style = 3)
        }
        for (i in 1:N){
            res[i, ] <- FUN(i, fit, args, call.freqs, seeds, prog = FALSE)
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
            cat("  |", rep(" ", 70), "|   0%, ", 0, " / ", N, "\n", sep = "",
                file = "prog.txt")
        }
        res <- t(parSapplyLB(cluster, 1:N, FUN, fit, args, call.freqs, seeds, prog))
        stopCluster(cluster)
        if (prog){
            unlink("prog.txt")
        }
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
