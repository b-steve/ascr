## Code to simulate the standard error correction procedure.

## Generate 10*n calls from n frogs, analyse assuming independence,
## then carry out correction procedure.

## Letting R know where everything is.
admbsecr.dir <- "~/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Results/secorrect/1", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", sep = sep)

## Get required libraries.
library(parallel)
library(secr)
library(admbsecr)

## Loading trap positions.
setwd(dat.dir)
mics <- read.csv(file = "array1a-top.csv")
micnames <- 1:dim(mics)[1]
mics <- cbind(micnames, mics)
traps <- read.traps(data = mics, detector = "signal")

setwd(work.dir)
## Setup for simulations.
seed <- 15338
set.seed(seed)
nsims <- 500
nboots <- 1000
seeds <- sample(1:1000000, size = nsims)
buffer <- 35

## True parameter values.
D <- 1750
shape <- 2.39
scale <- 5.35
pars <- c(D = D, shape = shape, scale = scale)
detfn <- "th"
set.seed(seed)
## Calls per frog
cpf <- 10
bounds <- list(D = c(0, 20000), sigma <- c(0, 30))
scalefactors <- c(shape = 1000, scale = 1000)

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
cat("start: ", date(), "\n", file = "prog.txt", sep = "")

FUN <- function(i, traps, calls, mask, pars, detfn, bounds, scalefactors, nboots, seeds){
  set.seed(seeds[i])
  workdir <- getwd()
  dirname <- paste("fit", i, sep = ".")
  system(paste("mkdir", dirname, sep = " "))
  setwd(dirname)
  capt <- sim.capt(traps = traps, calls = calls, mask = mask, pars = pars, detfn = "th")
  fit <- try.admbsecr(sv = pars, capt = capt, traps = traps, mask = mask,
                      detfn = detfn, bounds = bounds, scalefactors = scalefactors)
  if (class(fit)[1] == "try-error"){
    fit.sec <- "error"
  } else {
    fit.sec <- se.correct(fit, size = nboots, calls = calls)
  }
  out <- list(fit = fit, fit.sec = fit.sec)
  setwd(workdir)
  system(paste("rm -rf", dirname, sep = " "))
  filename <- paste("fits/fits", i, "RData", sep = ".")
  save(out, file = filename)
  cat(i, ": ", date(), "\n", file = "prog.txt", append = TRUE)
  out
}

ncores <- getOption("cl.cores", detectCores())
myCluster <- makeCluster(ncores)
clusterEvalQ(myCluster, {
  require(admbsecr)
})
parallel.time <- system.time({
  res.parallel <- parLapply(myCluster, 1:nsims, FUN, traps = traps, calls = cpf,
                            mask = mask, pars = pars, detfn = detfn, bounds = bounds,
                            scalefactors = scalefactors, nboots = nboots,
                            seeds = seeds)
})
stopCluster(myCluster)
save.image(file = "allres.RData")
