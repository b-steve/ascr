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
seed <- 9578
D <- 5270
g0 <- 0.99999
sigma <- 5.60
sigmatoa <- 0.002
ssb0 <- 5.16
ssb1 <- -0.02
sigmass <- 6.20
cutoff <- 150
truepars <- c(D = D, g0 = g0, sigma = sigma, sigmatoa = sigmatoa,
              ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)

set.seed(seed)
## Inverse of speed of sound (in ms per metre).
invsspd <- 1000/330

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(traps), as.matrix(mask))
simprobs <- NULL
toaprobs <- NULL
ssprobs <- NULL
jointprobs <- NULL
simpleres <- matrix(0, nrow = nsims, ncol = 3)
toares <- matrix(0, nrow = nsims, ncol = 4)
ssres <- matrix(0, nrow = nsims, ncol = 4)
jointres <- matrix(0, nrow = nsims, ncol = 5)
colnames(simpleres) <- c("D", "g0", "sigma")
colnames(toares) <- c("D", "g0", "sigma", "sigmatoa")
colnames(ssres) <- c("D", "ssb0", "ssb1", "sigmass")
colnames(jointres) <- c("D", "sigmatoa", "ssb0", "ssb1", "sigmass")

## Carrying out simulation.
for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else {
    print(c(i, date()))
  }
  ## Simulating data and setting things up for analysis.
  popn <- sim.popn(D = D, core = traps, buffer = buffer)
  capthist.ss <- sim.capthist.ss(traps, popn, detectpars)
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
  ## Generating times of arrival.
  ## Time of call itself doesn't provide any information, this comes
  ## solely from the *differences* in arrival times between traps. We
  ## can assume that animal with ID = t calls at time t without loss
  ## of generality.
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
  ## Straightforward SECR model using admbsecr().
  simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                            sv = truepars[1:3], admbwd = admb.dir,
                            method = "simple", verbose = FALSE), silent = TRUE)
  if (class(simplefit) == "try-error"){
    simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                              sv = "auto", admbwd = admb.dir, method = "simple",
                              verbose = FALSE), silent = TRUE)
  }
  if (class(simplefit) == "try-error"){
    simplecoef <- NA
    simprobs <- c(simprobs, i)
  } else {
    simplecoef <- coef(simplefit)
  }
  ## SECR model using supplementary TOA data.
  ssqtoa <- apply(capthist.toa, 1, toa.ssq, dists = mask.dists)
  toafit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask,
                         sv = truepars[1:4], ssqtoa = ssqtoa, admbwd = admb.dir,
                         method = "toa", verbose = FALSE), silent = TRUE)
  if (class(toafit) == "try-error"){
    toafit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask,
                           sv = "auto", ssqtoa = ssqtoa, admbwd = admb.dir,
                           method = "toa", verbose = FALSE), silent = TRUE)
  }
  if (class(toafit) == "try-error"){
    toacoef <- NA
    toaprobs <- c(toaprobs, i)
  } else {
    toacoef <- coef(toafit)
  }
  ## SECR model using supplementary signal strength data.
  ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                        sv = truepars[c(1, 5:7)], cutoff = cutoff,
                        admbwd = admb.dir, method = "ss", verbose = FALSE),
               silent = FALSE)
  if (class(ssfit) == "try-error"){
    ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                          sv = "auto", cutoff = cutoff, admbwd = admb.dir,
                          method = "ss", verbose = FALSE), silent = FALSE)
  }
  if (class(ssfit) == "try-error"){
    sscoef <- NA
    ssprobs <- c(ssprobs, i)
  } else {
    sscoef <- coef(ssfit)
  }
  ## SECR model using joint TOA and signal strength data.
  jointfit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                           sv = truepars[c(1, 4:7)], cutoff = cutoff,
                           ssqtoa = ssqtoa, admbwd = admb.dir, method = "sstoa",
                           verbose = FALSE), silent = TRUE)
  if (class(jointfit) == "try-error"){
    jointfit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                             sv = "auto", cutoff = cutoff, ssqtoa = ssqtoa,
                             admbwd = admb.dir, method = "sstoa",
                             verbose = FALSE), silent = TRUE)
  }
  if (class(jointfit) == "try-error"){
    jointcoef <- NA
    jointprobs <- c(jointprobs, i)
  } else {
    jointcoef <- coef(jointfit)
  }
  simpleres[i, ] <- simplecoef
  toares[i, ] <- toacoef
  ssres[i, ] <- sscoef
  jointres[i, ] <- jointcoef
  if (i == nsims){
    print(c("end", date()))
  }
}

## To write the simulation results to a file.
write.table(simpleres, "~/admbsecr/Results/frogs/2/simpleres.txt", row.names = FALSE)
write.table(toares, "~/admbsecr/Results/frogs/2/toares.txt", row.names = FALSE)
write.table(ssres, "~/admbsecr/Results/frogs/2/ssres.txt", row.names = FALSE)
write.table(jointres, "~/admbsecr/Results/frogs/2/jointres.txt", row.names = FALSE)
