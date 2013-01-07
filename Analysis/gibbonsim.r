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

setwd(work.dir)
## Setup for simulations.
nsims <- 1
set.seed <- 7862
buffer <- 6000
mask.spacing <- 50
trap.spacing <- 500

## True parameter values.
D <- 0.0416
g0 <- 1
sigma <- 1250
kappa <- 70
truepars <- c(D = D, g0 = g0, sigma = sigma, kappa = kappa)
detectpars <- list(g0 = g0, sigma = sigma)

## Setting up mask and traps.
traps <- make.grid(nx = 3, ny = 1, spacing = trap.spacing, detector = "proximity")
ntraps <- nrow(traps)
mask <- make.mask(traps, spacing = mask.spacing, type = "trapbuffer", buffer = buffer)
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(traps), as.matrix(mask))
mask.angs <- angles(as.matrix(traps), as.matrix(mask))
simprobs <- NULL
angprobs <- NULL
simpleres <- matrix(0, nrow = nsims, ncol = 3)
angres <- matrix(0, nrow = nsims, ncol = 4)
colnames(simpleres) <- c("D", "g0", "sigma")
colnames(angres) <- c("D", "g0", "sigma", "kappa")

## Carrying out simulation.
for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else {
    print(c(i, date()))
  }
  ## Simulating data and setting things up for analysis.
  popn <- sim.popn(D = D, core = traps, buffer = buffer)
  capthist <- sim.capthist(traps, popn, detectfn = 0, detectpar = detectpars,
                           noccasions = 1, renumber = FALSE)
  n <- nrow(capthist)
  ndets <- sum(capthist)
  cue.ids <- unique(as.numeric(rownames(capthist)))
  detections <- popn[cue.ids, ]
  radians <- t(angles(as.matrix(traps), as.matrix(detections)))
  radians <- array(radians, c(dim(radians), 1))
  errors <- array(rvm(ntraps*n, mean = 0, k = kappa), dim(radians))
  radians <- (radians + errors) %% (2*pi)
  radians[radians == 0] <- 2*pi
  radhist <- capthist
  radhist[radhist == 1] <- radians[radhist == 1]
  ## Straightforward SECR model using admbsecr()
  simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                            #fix = list(g0 = 1),
                            sv = truepars[1:3], admbwd = admb.dir,
                            method = "simple", verbose = FALSE, autogen = FALSE),
                   silent = TRUE)
  if (class(simplefit)[1] == "try-error"){
    simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                              #fix = list(g0 = 1),
                              sv = "auto", admbwd = admb.dir,
                              method = "simple", verbose = FALSE, autogen = FALSE),
                     silent = TRUE)
  }
  if (class(simplefit)[1] == "try-error"){
    simplecoef <- NA
    simprobs <- c(simprobs, i)
  } else {
    simplecoef <- coef(simplefit)
  }
  ## SECR model using supplementary angle data
  angfit <- try(admbsecr(capt = radhist, traps = traps, mask = mask,
                         #fix = list(g0 = 1),
                         sv = truepars, admbwd = admb.dir,
                         method = "ang", verbose = FALSE, autogen = FALSE),
                silent = TRUE)
  if (class(angfit)[1] == "try-error"){
    angfit <- try(admbsecr(capt = radhist, traps = traps, mask = mask,
                           #fix = list(g0 = 1),
                           sv = "auto", admbwd = admb.dir,
                           method = "ang", verbose = FALSE, autogen = FALSE),
                  silent = TRUE)
  }
  if (class(angfit)[1] == "try-error"){
    angcoef <- NA
    angprobs <- c(angprobs, i)
  } else {
    angcoef <- coef(angfit)
  }
  simpleres[i, ] <- simplecoef
  angres[i, ] <- angcoef
  if (i == nsims){
    print(c("end", date()))
  }
}

##To write the simulation results to a file.
write.table(angres, "~/admbsecr/Results/gibbons/5/angres.txt", row.names = FALSE)
write.table(simpleres, "~/admbsecr/Results/gibbons/5/simpleres.txt", row.names = FALSE)

