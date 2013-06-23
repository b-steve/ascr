## Code to test parallelisation.

## Letting R know where everything is.
admbsecr.dir <- "~/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Analysis/tests", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", sep = sep)

## Get required libraries.
library(secr)

## Loading the admbsecr library.
library(parallel)
library(admbsecr)

## Loading trap positions.
setwd(dat.dir)
mics <- read.csv(file = "array1a-top.csv")
micnames <- 1:dim(mics)[1]
mics <- cbind(micnames, mics)
traps <- read.traps(data = mics, detector = "signal")
setwd(work.dir)

## Setup for simulations.
nsims <- 4
buffer <- 35

## True parameter values.
seed <- 3038
## D 10 times smaller than usual as each frog emits 10 sounds.
D <- 1750
shape <- 2.39
scale <- 5.35
pars <- c(D = D, shape = shape, scale = scale)
bounds <- NULL
set.seed(seed)
cpf <- 1
scalefactors <- c(shape = 1000, scale = 1000)

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")

FUN <- function(i, traps, calls, mask, pars, detfn, bounds, scalefactors){
  workdir <- getwd()
  dirname <- paste("fit", i, sep = ".")
  system(paste("mkdir", dirname, sep = " "))
  setwd(dirname)
  capt <- sim.capt(traps = traps, calls = calls, mask = mask, pars = pars, detfn = "th")
  fit <- try.admbsecr(sv = pars, capt = capt, traps = traps, mask = mask,
                      detfn = detfn, bounds = bounds, scalefactors = scalefactors)
  setwd(workdir)
  system(paste("rm -rf", dirname, sep = " "))
  fit
}

series.time <- system.time({
  res.series <- lapply(1:nsims, FUN, traps = traps, calls = cpf, mask = mask, pars = pars,
                       detfn = "th", bounds = bounds, scalefactors = scalefactors)
})

ncores <- getOption("cl.cores", detectCores())
myCluster <- makeCluster(ncores)
clusterEvalQ(myCluster, {
  require(admbsecr)
})

parallel.time <- system.time({
  res.parallel <- parLapply(myCluster, 1:nsims, FUN, traps = traps, calls = cpf,
                            mask = mask, pars = pars, detfn = "th", bounds = bounds,
                            scalefactors = scalefactors)
})
stopCluster(myCluster)
