## admbsecr() takes capture history and mask objects from the secr
## package and fits an SECR model using ADMB.
admbsecr <- function(capt, traps, mask, sv = c(2000, 0.9, 10, 5), ssqtoa = NULL, angs = NULL, admbwd = NULL, method = "simple"){
  require(R2admb)
  require(secr)
  ## Warnings for incurrect input.
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
  ## Extracting no. animals trapped (n) and traps (k) from capture history array.
  ## Only currently works with one capture session.
  n <- dim(capt)[1]
  k <- dim(capt)[3]
  ## Area covered by each mask location.
  A <- attr(mask, "area")
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
    ang <- angles(traps, mask)
    data <- list(n = n, ntraps = k, nmask = nm, A = A, angcapt = capt, ang = ang, dist = dist)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], kappa = sv[4])
    bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000), kappa = c(0, 100000))
  } else {
    stop('method must be either "simple" or "toa"')
  }
  ## Adding bounds to the g0 parameter.
  
  ## Fitting the model.
  fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = TRUE,
                 run.opts = run.control(checkdata = "write", checkparam = "write", clean = TRUE))
  setwd(currwd)
  fit
}
