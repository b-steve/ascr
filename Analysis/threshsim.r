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

## True parameter values.
seed <- 8572
D <- 1700
ssb0 <- 160
ssb1 <- -2.40
sigmass <- 9.0
sigmatoa <- 0.0018
cutoff <- 130
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)
## Sensible start values for different detection functions.
hn.start <- c(D = D, g0 = 0.99999, sigma = 5.60)
hr.start <- c(D = D, g0 = 0.99999, sigma = 5.60, z = 12)
th.start <- c(D = D, shape = -3.70, scale = -0.51)
logth.start <- c(D = D, shape1 = 2.29, shape2 = 3.38, scale = -0.36)
ss.start <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
logss.start <- c(D = D, ssb0 = log(ssb0), ssb1 = -0.1, sigmass = sigmass)
ssmrds.start <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
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

hnres <- matrix(0, nrow = nsims, ncol = 6)
hnfixres <- matrix(0, nrow = nsims, ncol = 5)
hrres <- matrix(0, nrow = nsims, ncol = 7)
hrfixres <- matrix(0, nrow = nsims, ncol = 6)
thres <- matrix(0, nrow = nsims, ncol = 6)
##logthres <- matrix(0, nrow = nsims, ncol = 7)
ssres <- matrix(0, nrow = nsims, ncol = 7)
logssres <- matrix(0, nrow = nsims, ncol = 7)
mrdsres <- matrix(0, nrow = nsims, ncol = 6)
ssmrdsres <- matrix(0, nrow = nsims, ncol = 7)
toahnres <- matrix(0, nrow = nsims, ncol = 7)
toahnfixres <- matrix(0, nrow = nsims, ncol = 6)
toahrres <- matrix(0, nrow = nsims, ncol = 8)
toahrfixres <- matrix(0, nrow = nsims, ncol = 7)
toathres <- matrix(0, nrow = nsims, ncol = 7)
sstoares <- matrix(0, nrow = nsims, ncol = 8)
logsstoares <- matrix(0, nrow = nsims, ncol = 8)

colnames(hnres) <- c("D", "g0", "sigma", "logLik", "AIC", "maxgrad")
colnames(hnfixres) <- c("D", "sigma", "logLik", "AIC", "maxgrad")
colnames(hrres) <- c("D", "g0", "sigma", "z", "logLik", "AIC", "maxgrad")
colnames(hrfixres) <- c("D", "sigma", "z", "logLik", "AIC", "maxgrad")
colnames(thres) <- c("D", "shape", "scale", "logLik", "AIC", "maxgrad")
##colnames(logthres) <- c("D", "shape1", "shape2", "scale", "logLik", "AIC", "maxgrad")
colnames(ssres) <- c("D", "ssb0", "ssb1", "sigmass", "logLik", "AIC", "maxgrad")
colnames(logssres) <- c("D", "ssb0", "ssb1", "sigmass", "logLik", "AIC", "maxgrad")
colnames(mrdsres) <- c("D", "shape", "scale", "logLik", "AIC", "maxgrad")
colnames(ssmrdsres) <- c("D", "ssb0", "ssb1", "sigmass", "logLik", "AIC", "maxgrad")
colnames(toahnres) <- c("D", "g0", "sigma", "sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toahnfixres) <- c("D", "sigma", "sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toahrres) <-  c("D", "g0", "sigma", "z", "sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toahrfixres) <- c("D", "sigma", "z", "sigmatoa", "logLik", "AIC", "maxgrad")
colnames(toathres) <- c("D", "shape", "scale", "sigmatoa", "logLik", "AIC", "maxgrad")
colnames(sstoares) <- c("D", "ssb0", "ssb1", "sigmass", "sigmatoa", "logLik", "AIC", "maxgrad")
colnames(logsstoares) <- c("D", "ssb0", "ssb1", "sigmass", "sigmatoa", "logLik", "AIC", "maxgrad")
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
  hnfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                            sv = hn.start, method = "simple", detfn = "hn")
                   , silent = TRUE)
  if (class(hnfit)[1] == "try-error"){
    hnfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                              sv = "auto", method = "simple", detfn = "hn")
                     , silent = TRUE)
  }
  if (!class(hnfit)[1] == "try-error"){
    hncoef <- c(coef(hnfit), logLik(hnfit), AIC(hnfit), hnfit$maxgrad)
  } else hncoef <- NA
  ## TOA with half normal detection function.
  toahnfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask, bounds = bounds,
                           sv = c(hn.start, sigmatoa = sigmatoa), method = "toa", detfn = "hn")
                  , silent = TRUE)
  if (class(toahnfit)[1] == "try-error"){
    toahnfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask, bounds = bounds,
                             sv = "auto", method = "toa", detfn = "hn")
                    , silent = TRUE)
  }
  if (!class(hnfit)[1] == "try-error"){
    toahncoef <- c(coef(toahnfit), logLik(toahnfit), AIC(toahnfit), toahnfit$maxgrad)
  } else toahncoef <- NA
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
  if (!class(hnfixfit)[1] == "try-error"){
    hnfixcoef <- c(coef(hnfixfit), logLik(hnfixfit), AIC(hnfixfit), hnfixfit$maxgrad)
  } else hnfixcoef <- NA
  ## TOA with half normal detection function with fixed g0.
  toahnfixfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask, bounds = bounds,
                              sv = c(hn.start[-2], sigmatoa = sigmatoa), fix = list(g0 = 1),
                              method = "toa", detfn = "hn"), silent = TRUE)
  if (class(toahnfixfit)[1] == "try-error"){
    toahnfixfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask, bounds = bounds,
                             sv = "auto", fix = list(g0 = 1), method = "toa", detfn = "hn")
                       , silent = TRUE)
  }
  if (!class(hnfit)[1] == "try-error"){
    toahnfixcoef <- c(coef(toahnfixfit), logLik(toahnfixfit), AIC(toahnfixfit), toahnfixfit$maxgrad)
  } else toahnfixcoef <- NA
  ## Hazard rate detection function.
  hrfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                        bounds = bounds, sv = hr.start, method = "simple", detfn = "hr")
               , silent = TRUE)
  if (!class(hrfit)[1] == "try-error"){
    hrcoef <- c(coef(hrfit), logLik(hrfit), AIC(hrfit), hrfit$maxgrad)
  } else hrcoef <- NA
  ## TOA with hazard rate detection function.
  toahrfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask,
                           bounds = bounds, sv = c(hr.start, sigmatoa = sigmatoa),
                           method = "toa", detfn = "hr")
                  , silent = TRUE)
  if (!class(hrfit)[1] == "try-error"){
    toahrcoef <- c(coef(toahrfit), logLik(toahrfit), AIC(toahrfit), toahrfit$maxgrad)
  } else toahrcoef <- NA
  ## Hazard rate detection function with fixed g0.
  hrfixfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                           sv = hr.start[-2], fix = list(g0 = 1),
                           bounds = bounds, method = "simple", detfn = "hr")
                  , silent = TRUE)
  if (!class(hrfixfit)[1] == "try-error"){
    hrfixcoef <- c(coef(hrfixfit), logLik(hrfixfit), AIC(hrfixfit), hrfixfit$maxgrad)
  } else hrfixcoef <- NA
  ## TOA with hazard rate detection function with fixed g0.
  toahrfixfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask,
                              sv = c(hr.start[-2], sigmatoa = sigmatoa), fix = list(g0 = 1),
                              bounds = bounds, method = "toa", detfn = "hr")
                     , silent = TRUE)
  if (!class(hrfixfit)[1] == "try-error"){
    toahrfixcoef <- c(coef(toahrfixfit), logLik(toahrfixfit), AIC(toahrfixfit), toahrfixfit$maxgrad)
  } else toahrfixcoef <- NA
  ## Threshold detection function.
  thfit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, bounds = bounds,
                        sv = th.start, method = "simple", detfn = "th")
               , silent = TRUE)
  if (!class(thfit)[1] == "try-error"){
    thcoef <- c(coef(thfit), logLik(thfit), AIC(thfit), thfit$maxgrad)
  } else thcoef <- NA
  ## TOA with threshold detection function.
  toathfit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask, bounds = bounds,
                           sv = c(th.start, sigmatoa = sigmatoa), method = "toa", detfn = "th")
                  , silent = TRUE)
  if (!class(thfit)[1] == "try-error"){
    toathcoef <- c(coef(toathfit), logLik(toathfit), AIC(toathfit), toathfit$maxgrad)
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
  ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                        sv = ss.start, cutoff = cutoff, method = "ss",
                        bounds = bounds, detfn = "identity"), silent = TRUE)
  if (class(ssfit)[1] == "try-error"){
    ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                          sv = "auto", cutoff = cutoff, method = "ss",
                          bounds = bounds, detfn = "identity"), silent = TRUE)
  }
  if (!class(ssfit)[1] == "try-error"){
    sscoef <- c(coef(ssfit), logLik(ssfit), AIC(ssfit), ssfit$maxgrad)
  } else sscoef <- NA
  ## Joint TOA/SS analysis.
  sstoafit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                           sv = c(ss.start, sigmatoa = sigmatoa), cutoff = cutoff,
                           method = "sstoa", bounds = bounds, detfn = "identity")
                  , silent = TRUE)
  if (class(sstoafit)[1] == "try-error"){
    sstoafit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                             sv = "auto", cutoff = cutoff, method = "sstoa",
                             bounds = bounds, detfn = "identity"), silent = TRUE)
  }
  if (!class(sstoafit)[1] == "try-error"){
    sstoacoef <- c(coef(sstoafit), logLik(sstoafit), AIC(sstoafit), sstoafit$maxgrad)
  } else sstoacoef <- NA
  ## Signal strength detection function with log-link.
  logssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                           sv = logss.start, cutoff = cutoff, method = "ss",
                           detfn = "log", bounds = bounds),
                  silent = TRUE)
  if (class(logssfit)[1] == "try-error"){
    logssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                             bounds = bounds, sv = "auto", cutoff = cutoff, method = "ss",
                             detfn = "log"), silent = TRUE)
  }
  if (!class(logssfit)[1] == "try-error"){
    logsscoef <- c(coef(logssfit), logLik(logssfit), AIC(logssfit), logssfit$maxgrad)
  } else logsscoef <- NA
  ## Joint TOA/SS analysis with log-link.
  logsstoafit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                              sv = c(logss.start, sigmatoa = sigmatoa), cutoff = cutoff,
                              method = "sstoa", detfn = "log", bounds = bounds),
                     silent = TRUE)
  if (class(logsstoafit)[1] == "try-error"){
    logsstoafit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                                bounds = bounds, sv = "auto", cutoff = cutoff, method = "sstoa",
                                detfn = "log"), silent = TRUE)
  }
  if (!class(logsstoafit)[1] == "try-error"){
    logsstoacoef <- c(coef(logsstoafit), logLik(logsstoafit), AIC(logsstoafit), logsstoafit$maxgrad)
  } else logsstoacoef <- NA
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
  if (!class(mrdsfit)[1] == "try-error"){
    mrdscoef <- c(coef(mrdsfit), logLik(mrdsfit), AIC(mrdsfit), mrdsfit$maxgrad)
  } else mrdscoef <- NA
  ## MRDS with signal strength detection function.
  capthist.ssmrds <- array(0, dim = c(n, 1, 6, 2))
  capthist.ssmrds[, , , 1] <- capthist.ss
  capthist.ssmrds[, , , 2] <- dists
  ssmrdsfit <- try(admbsecr(capthist.ssmrds, traps = traps, mask = mask,
                            bounds = list(ssb0 = c(cutoff, 250)),
                            sv = ssmrds.start, cutoff = cutoff, method = "ssmrds",
                            trace = TRUE), silent = TRUE)
  if (!class(ssmrdsfit)[1] == "try-error"){
    ssmrdscoef <- c(coef(ssmrdsfit), logLik(ssmrdsfit), AIC(ssmrdsfit), ssmrdsfit$maxgrad)
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

