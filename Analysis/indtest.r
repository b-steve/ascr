## Code to test the violation of the independence assumpton.

## Generate 10*n calls from n frogs, analyse assuming independence and
## see what happens. Might also do some stuff on getting animal
## density from call density.

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
nsims <- 100
buffer <- 35

## True parameter values.
seed <- 3038
## D 10 times smaller than usual as each frog emits 10 sounds.
D <- 170
ssb0 <- 160
ssb1 <- 2.40
sigmass <- 9.0
sigmatoa <- 0.0018
cutoff <- 130
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)
sv <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass, sigmatoa = sigmatoa)
bounds <- NULL
set.seed(seed)
## Inverse of speed of sound (in ms per metre).
invsspd <- 1000/330
## Calls per frog
cpf <- 10

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")

## Function for stable model fitting.
try.admbsecr <- function(sv = "auto", ...){
  res <- try(admbsecr(sv = sv, ...), silent = TRUE)
  if (class(res)[1] == "try-error"){
    res <- try(admbsecr(sv = "auto", ...), silent = TRUE)
  }
  res
}

res <- matrix(0, nrow = nsims, ncol = 8)
colnames(res) <- c("D", "ssb0", "ssb1", "sigmass", "sigmatoa", "se.D",
                   "logLik", "maxgrad")

for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else {
    print(c(i, date()))
  }
  popn.i <- as.matrix(sim.popn(D = D, core = traps, buffer = buffer))
  popn <- matrix(0, nrow = cpf*nrow(popn.i), ncol = 2)
  for (j in 1:nrow(popn.i)){
    popn[(10*j - 9):(10*j), 1] <- popn.i[j, 1]
    popn[(10*j - 9):(10*j), 2] <- popn.i[j, 2]
  }
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
  fit <- try.admbsecr(sv = sv, capt = capthist.joint, traps = traps, mask = mask,
                      cutoff = cutoff, bounds = bounds, method = "sstoa",
                      detfn = "identity")
  if (class(fit)[1] != "try-error"){
    res[i, ] <- c(coef(fit), stdEr(fit)[1], logLik(fit), fit$maxgrad)
  } else{
    res[i, ] <- NA
  } 
  if (i == nsims){
    print(c("end", date()))
  }
}

write.table(res, file = "~/admbsecr/Results/ind/1/res.txt")
