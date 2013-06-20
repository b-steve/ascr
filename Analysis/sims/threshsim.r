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
library(R2admb)
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

## True parameter values.
seed <- 8572
D <- 1700
ssb0 <- 160
ssb1 <- 2.40
sigmass <- 9.0
sigmatoa <- 0.0018
cutoff <- 130
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)
## Sensible start values for different detection functions.
hn.start <- c(D = D, g0 = 0.99999, sigma = 8.20)
hr.start <- c(D = D, g0 = 0.99999, sigma = 9.50, z = 6.10)
th.start <- c(D = D, shape = 2.00, scale = 5.9)
logth.start <- c(D = D, shape1 = 2.29, shape2 = 3.38, scale = -0.36)
ss.start <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
logss.start <- c(D = D, ssb0 = log(ssb0), ssb1 = 0.05, sigmass = sigmass)
ssmrds.start <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
## Bounding D.
bounds <- list(D = c(0, 15000))
hr.bounds <- list(D = c(0, 15000), z = c(0, 20), sigma = c(1, 50))
toahr.bounds <- list(D = c(0, 15000), z = c(0, 20), sigmatoa = c(0, 0.1))
logsstoa.bounds <- list(D = c(0, 15000), ssb0 = c(0, 7), ssb1 = c(0, 5),
                        sigmatoa = c(0, 0.1), sigmass = c(0, 50))
set.seed(seed)
## Inverse of speed of sound (in ms per metre).
invsspd <- 1000/330

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(traps), as.matrix(mask))

## Function for stable model fitting.
try.admbsecr <- function(sv = "auto", ...){
  res <- try(admbsecr(sv = sv, ...), silent = TRUE)
  if (class(res)[1] == "try-error"){
    res <- try(admbsecr(sv = "auto", ...), silent = TRUE)
  }
  res
}

hnres <- matrix(0, nrow = nsims, ncol = 9)
hnfixres <- matrix(0, nrow = nsims, ncol = 7)
hrres <- matrix(0, nrow = nsims, ncol = 11)
hrfixres <- matrix(0, nrow = nsims, ncol = 9)
thres <- matrix(0, nrow = nsims, ncol = 9)
##logthres <- matrix(0, nrow = nsims, ncol = 11)
ssres <- matrix(0, nrow = nsims, ncol = 11)
logssres <- matrix(0, nrow = nsims, ncol = 11)
mrdsres <- matrix(0, nrow = nsims, ncol = 9)
ssmrdsres <- matrix(0, nrow = nsims, ncol = 11)
toahnres <- matrix(0, nrow = nsims, ncol = 11)
toahnfixres <- matrix(0, nrow = nsims, ncol = 9)
toahrres <- matrix(0, nrow = nsims, ncol = 13)
toahrfixres <- matrix(0, nrow = nsims, ncol = 11)
toathres <- matrix(0, nrow = nsims, ncol = 11)
sstoares <- matrix(0, nrow = nsims, ncol = 13)
logsstoares <- matrix(0, nrow = nsims, ncol = 13)

colnames(hnres) <- c("D", "g0", "sigma", "se.D", "se.g0", "se.sigma", "logLik", "AIC", "maxgrad")
colnames(hnfixres) <- c("D", "sigma", "se.D", "se.sigma", "logLik", "AIC", "maxgrad")
colnames(hrres) <- c("D", "g0", "sigma", "z", "se.D", "se.g0", "se.sigma", "se.z", "logLik", "AIC", "maxgrad")
colnames(hrfixres) <- c("D", "sigma", "z", "se.D", "se.sigma", "se.z", "logLik", "AIC", "maxgrad")
colnames(thres) <- c("D", "shape", "scale", "se.D", "se.shape", "se.scale", "logLik", "AIC", "maxgrad")
##colnames(logthres) <- c("D", "shape1", "shape2", "scale", "se.D", "se.shape1", "se.shape2", "se.scale", "logLik", "AIC", "maxgrad")
colnames(ssres) <- c("D", "ssb0", "ssb1", "sigmass", "se.D", "se.ssb0", "se.ssb1", "se.sigmass", "logLik", "AIC", "maxgrad")
colnames(logssres) <- c("D", "ssb0", "ssb1", "sigmass", "se.D", "se.ssb0", "se.ssb1", "se.sigmass", "logLik", "AIC", "maxgrad")
colnames(mrdsres) <- c("D", "shape", "scale", "se.D", "se.shape", "se.scale", "logLik", "AIC", "maxgrad")
colnames(ssmrdsres) <- c("D", "ssb0", "ssb1", "sigmass", "se.D", "se.ssb0", "se.ssb1", "se.sigmass", "logLik", "AIC", "maxgrad")
colnames(toahnres) <- c("D", "g0", "sigma", "sigmatoa", "se.D", "se.g0", "se.sigma", "se.sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toahnfixres) <- c("D", "sigma", "sigmatoa", "se.D", "se.sigma", "se.sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toahrres) <-  c("D", "g0", "sigma", "z", "sigmatoa", "se.D", "se.g0", "se.sigma", "se.z", "se.sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toahrfixres) <- c("D", "sigma", "z", "sigmatoa", "se.D", "se.sigma", "se.z", "se.sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toathres) <- c("D", "shape", "scale", "sigmatoa", "se.D", "se.shape", "se.scale", "se.sigmatoa", "logLik", "AIC", "maxgrad")
colnames(sstoares) <- c("D", "ssb0", "ssb1", "sigmass", "sigmatoa", "se.D", "se.ssb0", "se.ssb1", "se.sigmass", "se.sigmatoa", "logLik", "AIC", "maxgrad")
colnames(logsstoares) <- c("D", "ssb0", "ssb1", "sigmass", "sigmatoa", "se.D", "se.ssb0", "se.ssb1", "se.sigmass", "se.sigmatoa", "logLik", "AIC", "maxgrad")
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
  ## Generating TOA data (see frogsim.r)
  capthist.toa <- array(0, dim = dim(capthist))
  for (j in 1:n){
    for (k in 1:ntraps){
      if (capthist[j, 1, k] == 1){
        dist <- dists[j, k]
        meantoa <- cue.ids[j] + invsspd*dist/1000
        capthist.toa[j, 1, k] <- rnorm(1, meantoa, sigmatoa)
      } else {
        capthist.toa[j, 1, k] <- 0
      }
    }
  }
  capthist.joint <- array(c(capthist.ss, capthist.toa), dim = c(dim(capthist), 2))
  ## Half normal detection function.
  hnfit <- try.admbsecr(sv = hn.start, capt = capthist, traps = traps, mask = mask,
                        bounds = bounds, method = "simple", detfn = "hn")
  if (!class(hnfit)[1] == "try-error"){
    hncoef <- c(coef(hnfit), stdEr(hnfit), logLik(hnfit), AIC(hnfit), hnfit$maxgrad)
  } else hncoef <- NA
  ## TOA with half normal detection function.
  toahnfit <- try.admbsecr(sv = c(hn.start, sigmatoa = sigmatoa), capt = capthist.toa,
                           traps = traps, mask = mask, bounds = bounds, method = "toa",
                           detfn = "hn")
  if (!class(hnfit)[1] == "try-error"){
    toahncoef <- c(coef(toahnfit), stdEr(toahnfit), logLik(toahnfit), AIC(toahnfit), toahnfit$maxgrad)
  } else toahncoef <- NA
  ## Half normal detection function with fixed g0.
  hnfixfit <- try.admbsecr(sv = hn.start[-2], capt = capthist, traps = traps, mask = mask,
                           fix = list(g0 = 1), bounds = bounds, method = "simple", detfn = "hn")
  if (!class(hnfixfit)[1] == "try-error"){
    hnfixcoef <- c(coef(hnfixfit), stdEr(hnfixfit), logLik(hnfixfit), AIC(hnfixfit), hnfixfit$maxgrad)
  } else hnfixcoef <- NA
  ## TOA with half normal detection function with fixed g0.
  toahnfixfit <- try.admbsecr(sv = c(hn.start[-2], sigmatoa = sigmatoa), capt = capthist.toa,
                              traps = traps, mask = mask, bounds = bounds, fix = list(g0 = 1),
                              method = "toa", detfn = "hn")
  if (!class(hnfit)[1] == "try-error"){
    toahnfixcoef <- c(coef(toahnfixfit), stdEr(toahnfixfit), logLik(toahnfixfit), AIC(toahnfixfit), toahnfixfit$maxgrad)
  } else toahnfixcoef <- NA
  ## Hazard rate detection function.
  hrfit <- try.admbsecr(sv = hr.start, capt = capthist, traps = traps, mask = mask,
                        bounds = hr.bounds, method = "simple", detfn = "hr")
  if (!class(hrfit)[1] == "try-error"){
    hrcoef <- c(coef(hrfit), stdEr(hrfit), logLik(hrfit), AIC(hrfit), hrfit$maxgrad)
  } else hrcoef <- NA
  ## TOA with hazard rate detection function.
  toahrfit <- try.admbsecr(sv = c(hr.start, sigmatoa = sigmatoa), capt = capthist.toa,
                           traps = traps, mask = mask, bounds = hr.bounds,
                           method = "toa", detfn = "hr")
  if (!class(hrfit)[1] == "try-error"){
    toahrcoef <- c(coef(toahrfit), stdEr(toahrfit), logLik(toahrfit), AIC(toahrfit), toahrfit$maxgrad)
  } else toahrcoef <- NA
  ## Hazard rate detection function with fixed g0.
  hrfixfit <- try.admbsecr(sv = hr.start[-2], capt = capthist, traps = traps, mask = mask,
                           fix = list(g0 = 1), bounds = hr.bounds, method = "simple", detfn = "hr")
  if (!class(hrfixfit)[1] == "try-error"){
    hrfixcoef <- c(coef(hrfixfit), stdEr(hrfixfit), logLik(hrfixfit), AIC(hrfixfit), hrfixfit$maxgrad)
  } else hrfixcoef <- NA
  ## TOA with hazard rate detection function with fixed g0.
  toahrfixfit <- try.admbsecr(sv = c(hr.start[-2], sigmatoa = sigmatoa), capt = capthist.toa,
                              traps = traps, mask = mask, fix = list(g0 = 1),
                              bounds = toahr.bounds, method = "toa", detfn = "hr")
  if (!class(hrfixfit)[1] == "try-error"){
    toahrfixcoef <- c(coef(toahrfixfit), stdEr(toahrfixfit), logLik(toahrfixfit), AIC(toahrfixfit), toahrfixfit$maxgrad)
  } else toahrfixcoef <- NA
  ## Threshold detection function.
  thfit <- try.admbsecr(sv = th.start, capt = capthist, traps = traps, mask = mask, bounds = bounds,
                        method = "simple", detfn = "th")
  if (!class(thfit)[1] == "try-error"){
    thcoef <- c(coef(thfit), stdEr(thfit), logLik(thfit), AIC(thfit), thfit$maxgrad)
  } else thcoef <- NA
  ## TOA with threshold detection function.
  toathfit <- try.admbsecr(sv = c(th.start, sigmatoa = sigmatoa), capt = capthist.toa, traps = traps,
                           mask = mask, bounds = bounds, method = "toa", detfn = "th",
                           scalefactors = c(sigmatoa = 5000))
  if (!class(thfit)[1] == "try-error"){
    toathcoef <- c(coef(toathfit), stdEr(toathfit), logLik(toathfit), AIC(toathfit), toathfit$maxgrad)
  } else toathcoef <- NA
  ## Log-link threshold detection function.
  ## logthfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
  ##                          sv = logth.start, method = "simple", detfn = "logth",
  ##                          bounds = bounds)
  ##                 , silent = TRUE)
  ## if (!class(logthfit)[1] == "try-error"){
  ##   logthcoef <- c(coef(logthfit), logLik(logthfit), AIC(logthfit), logthfit$maxgrad)
  ## }
  ## Signal strength detection function.
  ssfit <- try.admbsecr(sv = ss.start, capt = capthist.ss, traps = traps, mask = mask,
                        cutoff = cutoff, method = "ss", bounds = bounds, detfn = "identity")
  if (!class(ssfit)[1] == "try-error"){
    sscoef <- c(coef(ssfit), stdEr(ssfit), logLik(ssfit), AIC(ssfit), ssfit$maxgrad)
  } else sscoef <- NA
  ## Joint TOA/SS analysis.
  sstoafit <- try.admbsecr(sv = c(ss.start, sigmatoa = sigmatoa), capt = capthist.joint,
                           traps = traps, mask = mask, cutoff = cutoff,
                           method = "sstoa", bounds = bounds, detfn = "identity")
  if (!class(sstoafit)[1] == "try-error"){
    sstoacoef <- c(coef(sstoafit), stdEr(sstoafit), logLik(sstoafit), AIC(sstoafit), sstoafit$maxgrad)
  } else sstoacoef <- NA
  ## Signal strength detection function with log-link.
  logssfit <- try.admbsecr(sv = logss.start, capt = capthist.ss, traps = traps, mask = mask,
                           cutoff = cutoff, method = "ss", detfn = "log", bounds = bounds)
  if (!class(logssfit)[1] == "try-error"){
    logsscoef <- c(coef(logssfit), stdEr(logssfit), logLik(logssfit), AIC(logssfit), logssfit$maxgrad)
  } else logsscoef <- NA
  ## Joint TOA/SS analysis with log-link.
  logsstoafit <- try.admbsecr(sv = c(logss.start, sigmatoa = sigmatoa), capt = capthist.joint,
                              traps = traps, mask = mask, cutoff = cutoff, method = "sstoa",
                              detfn = "log", bounds = logsstoa.bounds)
  if (!class(logsstoafit)[1] == "try-error"){
    logsstoacoef <- c(coef(logsstoafit), stdEr(logsstoafit), logLik(logsstoafit), AIC(logsstoafit), logsstoafit$maxgrad)
  } else logsstoacoef <- NA
  ## MRDS with threshold detection function.
  capthist.mrds <- array(0, dim = c(n, 1, 6, 2))
  capthist.mrds[, , , 1] <- capthist
  capthist.mrds[, , , 2] <- dists
  mrdsfit <- try.admbsecr(sv = coef(thfit)[1:3], capthist.mrds, traps = traps, mask = mask,
                          bounds = list(scale = c(1, 15)), method = "mrds", detfn = "th")
  if (!class(mrdsfit)[1] == "try-error"){
    mrdscoef <- c(coef(mrdsfit), stdEr(mrdsfit), logLik(mrdsfit), AIC(mrdsfit), mrdsfit$maxgrad)
  } else mrdscoef <- NA
  ## MRDS with signal strength detection function.
  capthist.ssmrds <- array(0, dim = c(n, 1, 6, 2))
  capthist.ssmrds[, , , 1] <- capthist.ss
  capthist.ssmrds[, , , 2] <- dists
  ssmrdsfit <- try.admbsecr(sv = ssmrds.start, capthist.ssmrds, traps = traps, mask = mask,
                            bounds = list(ssb0 = c(cutoff, 250)), cutoff = cutoff, method = "ssmrds")
  if (!class(ssmrdsfit)[1] == "try-error"){
    ssmrdscoef <- c(coef(ssmrdsfit), stdEr(ssmrdsfit), logLik(ssmrdsfit), AIC(ssmrdsfit), ssmrdsfit$maxgrad)
  } else ssmrdscoef <- NA
  hnres[i, ] <- hncoef
  toahnres[i, ] <- toahncoef
  hnfixres[i, ] <- hnfixcoef
  toahnfixres[i, ] <- toahnfixcoef
  hrres[i, ] <- hrcoef
  toahrres[i, ] <- toahrcoef
  hrfixres[i, ] <- hrfixcoef
  toahrfixres[i, ] <- toahrfixcoef
  thres[i, ] <- thcoef
  toathres[i, ] <- toathcoef
  ##logthres[i, ] <- logthcoef
  ssres[i, ] <- sscoef
  sstoares[i, ] <- sstoacoef
  logssres[i, ] <- logsscoef
  logsstoares[i, ] <- logsstoacoef
  mrdsres[i, ] <- mrdscoef
  ssmrdsres[i, ] <- ssmrdscoef
  if (i == nsims){
    print(c("end", date()))
  }
}
## To write the simulation results to a file.
write.table(hnres, "~/admbsecr/Results/thresh/2/hnres.txt", row.names = FALSE)
write.table(hnfixres, "~/admbsecr/Results/thresh/2/hnfixres.txt", row.names = FALSE)
write.table(hrres, "~/admbsecr/Results/thresh/2/hrres.txt", row.names = FALSE)
write.table(hrfixres, "~/admbsecr/Results/thresh/2/hrfixres.txt", row.names = FALSE)
write.table(thres, "~/admbsecr/Results/thresh/2/thres.txt", row.names = FALSE)
##write.table(logthres, "~/admbsecr/Results/thresh/1/logthres.txt", row.names = FALSE)
write.table(ssres, "~/admbsecr/Results/thresh/2/ssres.txt", row.names = FALSE)
write.table(logssres, "~/admbsecr/Results/thresh/2/logssres.txt", row.names = FALSE)
write.table(mrdsres, "~/admbsecr/Results/thresh/2/mrdsres.txt", row.names = FALSE)
write.table(ssmrdsres, "~/admbsecr/Results/thresh/2/ssmrdsres.txt", row.names = FALSE)
write.table(toahnres, "~/admbsecr/Results/thresh/2/toahnres.txt", row.names = FALSE)
write.table(toahnfixres, "~/admbsecr/Results/thresh/2/toahnfixres.txt", row.names = FALSE)
write.table(toahrres, "~/admbsecr/Results/thresh/2/toahrres.txt", row.names = FALSE)
write.table(toahrfixres, "~/admbsecr/Results/thresh/2/toahrfixres.txt", row.names = FALSE)
write.table(toathres, "~/admbsecr/Results/thresh/2/toathres.txt", row.names = FALSE)
write.table(sstoares, "~/admbsecr/Results/thresh/2/sstoares.txt", row.names = FALSE)
write.table(logsstoares, "~/admbsecr/Results/thresh/2/logsstoares.txt", row.names = FALSE)
