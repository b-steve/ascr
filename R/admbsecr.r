## admbsecr() takes capture history and mask objects from the secr
## package and fits an SECR model using ADMB.
admbsecr <- function(capt, traps, mask, sv = "auto", ssqtoa = NULL,
                     angs = NULL, admbwd = NULL, method = "simple",
                     profpars = NULL){
    require(R2admb)
    require(secr)
    ## Warnings for incorrect input.
    if (length(method) != 1){
        stop("method must be of length 1")
    }
    if (method == "simple" & any(capt != 1 & capt != 0)){
        stop('capt must be binary when using the "simple" method')
    }
    prefix <- paste(method, "secr", sep = "")
    currwd <- getwd()
    ## Moving to ADMB working directory.
    if (!is.null(admbwd)){
        setwd(admbwd)
    }
    ## If NAs are present in capture history object, change to zeros.
    capt[is.na(capt)] <- 0
    ## Extracting no. animals trapped (n) and traps (k) from capture history array.
    ## Only currently works with one capture session.
    n <- dim(capt)[1]
    k <- dim(capt)[3]
    ## Area covered by each mask location.
    A <- attr(mask, "area")
    ## Setting sensible start values if elements of sv are "auto".
    if (length(sv) == 1 & sv[1] == "auto"){
        npars <- 3 + sum(method != "simple")
        sv <- rep("auto", npars)
        if (method == "ss"){
            sv <- c(sv, rep("auto", 2))
        }
    }
    if (any(sv == "auto")){
        ## Give sv vector names if it doesn't have them.
        if (is.null(names(sv))){
            names(sv) <- c("D", "g0", "sigma", "sigmatoa"[method == "toa"],
                           "kappa"[method == "ang"], "ssb0"[method == "ss"],
                           "ssb1"[method == "ss"], "sigmass"[method == "ss"])
        } else {
            ## Reordering sv vector if names are provided.
            sv <- sv[c("D", "g0", "sigma", "sigmatoa"[method == "toa"],
                       "kappa"[method == "ang"], "ssb0"[method == "ss"],
                       "ssb1"[method == "ss"], "sigmass"[method = "ss"])]
        }
        bincapt <- capt
        bincapt[capt > 0] <- 1
        autofuns <- list("D" = autoD, "g0" = autog0, "sigma" = autosigma,
                         "sigmatoa" = autosigmatoa, "kappa" = autokappa,
                         "ssb0" = autossb0, "ssb1" = autossb1,
                         "sigmass" = autosigmass)
        ## Replacing "auto" elements of sv vector.
        for (i in rev(which(sv == "auto"))){
            sv[i] <- autofuns[[names(sv)[i]]](capt, bincapt, traps, mask, sv)
        }
        sv <- as.numeric(sv)
    }
    ## Removing attributes from capt and mask objects as do_admb cannot handle them.
    capt <- matrix(as.vector(capt), nrow = n, ncol = k)
    mask <- as.matrix(mask)
    ## No. of mask locations.
    nm <- nrow(mask)
    ## Distances between traps and mask locations.
    dist <- distances(traps, mask)
    ## Setting up parameters for do_admb.
    if (method == "simple"){
        data <- list(n = n, ntraps = k, nmask = nm, A = A, capt = capt, dist = dist)
        params <- list(D = sv[1], g0 = sv[2], sigma = sv[3])
        bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000))
    } else if (method == "toa"){
        if (is.null(ssqtoa)){
            ssqtoa <- apply(capt, 1, toa.ssq, dists = dist)
        }
        data <- list(n = n, ntraps = k, nmask = nm, A = A, toacapt = capt,
                     toassq = t(ssqtoa), dist = dist)
        params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], sigmatoa = sv[4])
        bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000), sigmatoa = c(0, 100000))
    } else if (method == "ang"){
        if (is.null(angs)){
            angs <- angles(traps, mask)
        }
        data <- list(n = n, ntraps = k, nmask = nm, A = A, angcapt = capt, ang = angs, dist = dist)
        params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], kappa = sv[4])
        bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000), kappa = c(0, 100000))
    } else if (method == "ss"){
        data <- list(n = n, ntraps = k, nmask = nm, A = A, sscapt = capt, dist = dist)
        params <- list(D = sv[1], g0 = sv[2], sigma = sv[3],
                       ssb0 = sv[4], ssb1 = sv[5], sigmass = sv[6])
        bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000),
                       sigmass = c(0, 100000))
    } else {
        stop('method must be either "simple", "toa", "ang" or "ss"')
    }
    ## Fitting the model.
    if (!is.null(profpars)){
        fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = TRUE,
                       profile = TRUE, profpars = profpars,
                       run.opts = run.control(checkdata = "write", checkparam = "write", clean = TRUE))
    } else {
        fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = TRUE,
                       run.opts = run.control(checkdata = "write", checkparam = "write", clean = TRUE))
    }
    setwd(currwd)
    fit
}
