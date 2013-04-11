## Letting R know where everything is.
admbsecr.dir <- "~/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Analysis", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", sep = sep)

## Get required libraries.
library(secr)

## Loading the admbsecr library.
setwd(admbsecr.dir)
library(devtools)
load_all(".")
##library(admbsecr)

## Loading trap positions.
setwd(dat.dir)
mics <- read.csv(file = "array1a-top.csv")
micnames <- 1:dim(mics)[1]
mics <- cbind(micnames, mics)
traps <- read.traps(data = mics, detector = "signal")

setwd(work.dir)
## Setup for simulations.
nsims <- 500
buffer <- 35
mask.spacing <- 50

## True parameter values.
seed <- 6846
D <- 5270
ssb0 <- 173
ssb1 <- -3.1
sigmass <- 6.20
cutoff <- 150
truepars <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)
## Sensible start values for different detection functions.
hn.start <- c(D = D, g0 = 0.99999, sigma = 5.60)
hr.start <- c(D = D, g0 = 0.99999, sigma = 5.60, z = 12)
th.start <- c(D = D, shape = -3.70, scale = -0.51)
logth.start <- c(D = D, shape1 = 2.29, shape2 = 3.38, scale = -0.36)
ss.start <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
logss.start <- c(D = D, ssb0 = log(ssb0), ssb1 = -0.1, sigmass = sigmass)
## Bounding D.
bounds <- list(D = c(0, 15000))
set.seed(seed)
## Inverse of speed of sound (in ms per metre).
invsspd <- 1000/330

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(traps), as.matrix(mask))

hnprobs <- NULL
hnfixprobs <- NULL
hrprobs <- NULL
hrfixprobs <- NULL
thprobs <- NULL
logthprobs <- NULL
ssprobs <- NULL
logssprobs <- NULL
mrdsprobs <- NULL
hnres <- matrix(0, nrow = nsims, ncol = 6)
hnfixres <- matrix(0, nrow = nsims, ncol = 5)
hrres <- matrix(0, nrow = nsims, ncol = 7)
hrfixres <- matrix(0, nrow = nsims, ncol = 6)
thres <- matrix(0, nrow = nsims, ncol = 6)
logthres <- matrix(0, nrow = nsims, ncol = 7)
ssres <- matrix(0, nrow = nsims, ncol = 7)
logssres <- matrix(0, nrow = nsims, ncol = 7)
mrdsres <- matrix(0, nrow = nsims, ncol = 6)
colnames(hnres) <- c("D", "g0", "sigma", "logLik", "AIC", "maxgrad")
colnames(hnfixres) <- c("D", "sigma", "logLik", "AIC", "maxgrad")
colnames(hrres) <- c("D", "g0", "sigma", "z", "logLik", "AIC", "maxgrad")
colnames(hrfixres) <- c("D", "sigma", "z", "logLik", "AIC", "maxgrad")
colnames(thres) <- c("D", "shape", "scale", "logLik", "AIC", "maxgrad")
colnames(logthres) <- c("D", "shape1", "shape2", "scale", "logLik", "AIC", "maxgrad")
colnames(ssres) <- c("D", "ssb0", "ssb1", "sigmass", "logLik", "AIC", "maxgrad")
colnames(logssres) <- c("D", "ssb0", "ssb1", "sigmass", "logLik", "AIC", "maxgrad")
colnames(mrdsres) <- c("D", "shape", "scale", "logLik", "AIC", "maxgrad")

## Carrying out simulation.
for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else {
    print(c(i, date()))
  }
  ## Simulating data and setting things up for analysis.
  popn <- sim.popn(D = D, core = traps, buffer = buffer)
  capthist.ss <- sim.capthist.ss(traps, popn, detectpars, log.link = FALSE)
  capthist <- capthist.ss
  capthist[capthist > 0] <- 1
  n <- nrow(capthist)
  ndets <- sum(capthist)
  ## IDs for detected animals.
  cue.ids <- unique(as.numeric(rownames(capthist)))
  ## Cartesian coordinates of detected animals.
  detections <- popn[cue.ids, ]
  ## Distances from detected animals to traps.
  dists <- t(distances(as.matrix(traps), as.matrix(detections)))
  ## Half normal detection function.
  hnfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                            sv = hn.start, method = "simple", detfn = "hn")
                   , silent = TRUE)
  if (class(hnfit)[1] == "try-error"){
    hnfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                              sv = "auto", method = "simple", detfn = "hn")
                     , silent = TRUE)
  }
  if (class(hnfit)[1] == "try-error"){
    hncoef <- NA
    hnprobs <- c(hnprobs, i)
  } else {
    hncoef <- c(coef(hnfit), logLik(hnfit), AIC(hnfit), hnfit$maxgrad)
  }
  ## Half normal detection function with fixed g0.
  hnfixfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                            sv = hn.start[-2], fix = list(g0 = 1), bounds = bounds,
                            method = "simple", detfn = "hn")
                   , silent = TRUE)
  if (class(hnfixfit)[1] == "try-error"){
    hnfixfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                              sv = "auto", fix = list(g0 = 1), bounds = bounds,
                              method = "simple", detfn = "hn")
                     , silent = TRUE)
  }
  if (class(hnfixfit)[1] == "try-error"){
    hnfixcoef <- NA
    hnfixprobs <- c(hnfixprobs, i)
  } else {
    hnfixcoef <- c(coef(hnfixfit), logLik(hnfixfit), AIC(hnfixfit), hnfixfit$maxgrad)
  }
  ## Hazard rate detection function.
  hrfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                            sv = hr.start, method = "simple", detfn = "hr")
                   , silent = TRUE)
  if (class(hrfit)[1] == "try-error"){
    hrcoef <- NA
    hrprobs <- c(hrprobs, i)
  } else {
    hrcoef <- c(coef(hrfit), logLik(hrfit), AIC(hrfit), hrfit$maxgrad)
  }
  ## Half normal detection function with fixed g0.
  hrfixfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                           sv = hr.start[-2], fix = list(g0 = 1),
                           bounds = bounds, method = "simple",
                           detfn = "hr")
                  , silent = TRUE)
  if (class(hrfixfit)[1] == "try-error"){
    hrfixcoef <- NA
    hrfixprobs <- c(hrfixprobs, i)
  } else {
    hrfixcoef <- c(coef(hrfixfit), logLik(hrfixfit), AIC(hrfixfit), hrfixfit$maxgrad)
  }
  ## Threshold detection function.
  thfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                        sv = th.start, method = "simple", detfn = "th")
               , silent = TRUE)
  if (class(thfit)[1] == "try-error"){
    thcoef <- NA
    thprobs <- c(thprobs, i)
  } else {
    thcoef <- c(coef(thfit), logLik(thfit), AIC(thfit), thfit$maxgrad)
  }
  ## Log-link threshold detection function.
  ## logthfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
  ##                          sv = logth.start, method = "simple", detfn = "logth",
  ##                          bounds = bounds)
  ##                 , silent = TRUE)
  ## if (class(logthfit)[1] == "try-error"){
  ##   logthcoef <- NA
  ##   logthprobs <- c(logthprobs, i)
  ## } else {
  ##   logthcoef <- c(coef(logthfit), logLik(logthfit), AIC(logthfit), logthfit$maxgrad)
  ## }
  ## Signal strength detection function.
  ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                        sv = ss.start, cutoff = cutoff, method = "ss",
                        bounds = bounds,
                        detfn = "identity"), silent = TRUE)
  if (class(ssfit)[1] == "try-error"){
    ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                          sv = "auto", cutoff = cutoff, method = "ss",
                          bounds = bounds,
                          detfn = "identity"), silent = TRUE)
  }
  if (class(ssfit)[1] == "try-error"){
    sscoef <- NA
    ssprobs <- c(ssprobs, i)
  } else {
    sscoef <- c(coef(ssfit), logLik(ssfit), AIC(ssfit), ssfit$maxgrad)
  }
  ## Signal strength detection function with log-link.
  logssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                           sv = logss.start, cutoff = cutoff, method = "ss",
                           detfn = "log", bounds = bounds), silent = TRUE)
  if (class(logssfit)[1] == "try-error"){
    logssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                             bounds = bounds, sv = "auto", cutoff = cutoff, method = "ss",
                             detfn = "log"), silent = TRUE)
  }
  if (class(logssfit)[1] == "try-error"){
    logsscoef <- NA
    logssprobs <- c(logssprobs, i)
  } else {
    logsscoef <- c(coef(logssfit), logLik(logssfit), AIC(logssfit), logssfit$maxgrad)
  }
  ## MRDS with threshold detection function.
  capthist.mrds <- array(0, dim = c(n, 1, 6, 2))
  capthist.mrds[, , , 1] <- capthist
  capthist.mrds[, , , 2] <- dists
  
  mrdsfit <- try(admbsecr(capthist.mrds, traps = traps, mask = mask,
                           sv = coef(thfit)[2:3],
                           fix = list(D = 5270), method = "mrds",
                           detfn = "th"), silent = TRUE)
  mrdsfit <- try(admbsecr(capthist.mrds, traps = traps, mask = mask,
                          sv = c(D = 5270, coef(mrdsfit)), bounds = bounds,
                          method = "mrds", detfn = "th"), silent = TRUE)
  if (class(mrdsfit)[1] == "try-error"){
    mrdscoef <- NA
    mrdsprobs <- c(mrdsprobs, i)
  } else {
    mrdscoef <- c(coef(mrdsfit), logLik(mrdsfit), AIC(mrdsfit), mrdsfit$maxgrad)
  }
  hnres[i, ] <- hncoef
  hnfixres[i, ] <- hnfixcoef
  hrres[i, ] <- hrcoef
  hrfixres[i, ] <- hrfixcoef
  thres[i, ] <- thcoef
  ##logthres[i, ] <- logthcoef
  ssres[i, ] <- sscoef
  logssres[i, ] <- logsscoef
  mrdsres[i, ] <- mrdscoef
  if (i == nsims){
    print(c("end", date()))
  }
}

## To write the simulation results to a file.
write.table(hnres, "~/admbsecr/Results/thresh/1/hnres.txt", row.names = FALSE)
write.table(hnfixres, "~/admbsecr/Results/thresh/1/hnfixres.txt", row.names = FALSE)
write.table(hrres, "~/admbsecr/Results/thresh/1/hrres.txt", row.names = FALSE)
write.table(hrfixres, "~/admbsecr/Results/thresh/1/hrfixres.txt", row.names = FALSE)
write.table(thres, "~/admbsecr/Results/thresh/1/thres.txt", row.names = FALSE)
##write.table(logthres, "~/admbsecr/Results/thresh/1/logthres.txt", row.names = FALSE)
write.table(ssres, "~/admbsecr/Results/thresh/1/ssres.txt", row.names = FALSE)
write.table(logssres, "~/admbsecr/Results/thresh/1/logssres.txt", row.names = FALSE)
write.table(mrdsres, "~/admbsecr/Results/thresh/1/mrdsres.txt", row.names = FALSE)
