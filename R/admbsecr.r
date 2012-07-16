## admbsecr() takes capture history and mask objects from the secr
## package and fits an SECR model using ADMB.
admbsecr <- function(capt, mask, sv = c(0, 0.5, 1), admbwd = NULL, prefix){
  require(R2admb)
  require(secr)
  currwd <- getwd()
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
  data <- list(n = n, ntraps = k, nmask = nm, A = A, capt = capt, dist = dist)
  params <- list(logD = sv[1], g0 = sv[2], logsigma = sv[3])
  bounds <- list(g0 = c(0, 1))
  ## Fitting the model.
  fit <- do_admb(prefix, data = data, params = params,
                 run.opts = run.control(checkdata = "write"))
  setwd(currwd)
  fit
}


