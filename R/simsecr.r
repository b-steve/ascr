## Simulates capt matrix suitable for admbsecr(), from either a model
## fit (using estimated parameters) or from provided values.

sim.secr <- function(fit){
  method <- fit$method
  detfn <- fit$detfn
  traps <- fit$traps
  mask <- fit$mask
  fitcoefs <- coef(fit)
  allparnames <- fit$parnames
  coefparnames <- names(fitcoefs)
  dataparnames <- allparnames[!allparnames %in% coefparnames]
  fixcoefs <- numeric(length(dataparnames))
  allcoefs <- c(fitcoefs, fixcoefs)
  names(fixcoefs) <- dataparnames
  for (i in dataparnames){
    fixcoefs[i] <- fit$data[[i]]
  }
  buffer <- 0
  range.x <- range(mask[, 1])
  range.y <- range(mask[, 2])
  core <- data.frame(x = range.x, y = range.y)
  popn <- sim.popn(D = allcoefs["D"], core = core, buffer = buffer)
  popn.dists <- distances(as.matrix(popn), as.matrix(traps))
  bincapt <- sim.bincapt(popn.dists, method, detfn, allcoefs)
  ## Unfinished.
}

sim.bincapt <- function(popn.dists, method, detfn, allcoefs){
  ## Unfinished.
}
