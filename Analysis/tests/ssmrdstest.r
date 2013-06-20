## Simulating some data and testing the ssmrds method.

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

## Simulating some data to test ssmrds method.
buffer <- 35
mask.spacing <- 50
seed <- 1314
D <- 5270
ssb0 <- 170
ssb1 <- -2.50
sigmass <- 7.00
cutoff <- 150
truepars <- c(D = D, ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass,
                   cutval = cutoff)
set.seed(seed)
## Inverse of speed of sound (in ms per metre).
invsspd <- 1000/330
## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(traps), as.matrix(mask))

popn <- sim.popn(D = D, core = traps, buffer = buffer)
capthist.ss <- sim.capthist.ss(traps, popn, detectpars, log.link = FALSE)
n <- nrow(capthist.ss)
k <- nrow(traps)
cue.ids <- unique(as.numeric(rownames(capthist.ss)))
detections <- popn[cue.ids, ]
dists <- t(distances(as.matrix(traps), as.matrix(detections)))
capthist.ssmrds <- array(0, dim = c(n, 1, k, 2))
capthist.ssmrds[, , , 1] <- capthist.ss
capthist.ssmrds[, , , 2] <- dists
ssmrds.fit <- admbsecr(capt = capthist.ssmrds, traps = traps, mask = mask,
                       sv = "auto", cutoff = 150, method = "ssmrds",
                       clean = TRUE)
