## admbsecr() takes capture history and mask objects from the secr
## package and fits an SECR model using ADMB.
admbsecr <- function(capt, mask, sv = c(0, 0.5, 1, 1), admbwd = NULL, method = "simple"){
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
    params <- list(logD = sv[1], g0 = sv[2], logsigma = sv[3])
  } else if (method == "toa"){
    ssqtoa <- apply(capt, 1, toa.ssq, dists = dist)
    data <- list(n = n, ntraps = k, nmask = nm, A = A, toacapt = capt, toassq = ssqtoa, dist = dist)
    params <- list(logD = sv[1], g0 = sv[2], logsigma = sv[3], logsigmatoa = sv[4])
  } else {
    stop('method must be either "simple" or "toa"')
  }
  ## Adding bounds to the g0 parameter.
  bounds <- list(g0 = c(0, 1))
  ## Fitting the model.
  fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = TRUE,
                 run.opts = run.control(checkdata = "write", checkparam = "write"))
  setwd(currwd)
  fit
}
