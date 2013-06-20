## Code to generate data where individuals emit calls of varying
## source strength. Analysis uses methods "ss" and "sstoa", which
## assume identical source strengths.

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

## Get required library.
library(secr)
library(CircStats)

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
ntraps <- nrow(traps)
setwd(dat.dir)
buffer <- 35
mask.spacing <- 50
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")

## True parameter values.
D <- 5270
ssb0 <- 173
ssb1 <- -3.1
sigmass <- 0
sigmas <- 5.5
sigmatoa <- 0.002
cutoff <- 150
truepars <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass, sigmas = sigmas)
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, sds = sigmas, cutval = cutoff)
invsspd <- 1000/330

## Simulating population
popn <- sim.popn(D = D, core = traps, buffer = buffer)
## Simulating capture history
capthist.ss <- sim.capthist.ss(traps, popn, detectpars, log.link = FALSE, re = FALSE)
capthist <- capthist.ss
capthist[capthist > 0] <- 1
n <- nrow(capthist)
ndets <- sum(capthist)
cue.ids <- unique(as.numeric(rownames(capthist)))
detections <- popn[cue.ids, ]
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

simple.fit <- admbsecr(capt = capthist, traps = traps, mask = mask, sv = "auto",
                       method = "simple", trace = FALSE)

ss.fit <- admbsecr(capt = capthist.ss, traps = traps, mask = mask, sv = "auto",
                   cutoff = 150, method = "ss", trace = FALSE)

toa.fit <- admbsecr(capt = capthist.toa, traps = traps, mask  = mask, sv = "auto",
                    method = "toa", trace = FALSE)

joint.fit <- admbsecr(capt = capthist.joint, traps = traps, mask = mask, sv = "auto",
                      bounds = list(ssb0 = c(150, 200)),
                      method = "sstoa", cutoff = 150, trace = TRUE)
