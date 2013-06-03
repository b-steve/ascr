## Testing the angdist method, i.e., both angle and distance data
## provided. Using parameters motivated by the gibbon data.

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
##library(admbsecr) # Use this if not loading with devtools.

setwd(work.dir)
## Setup for simulations.
nsims <- 500
set.seed(3788)
buffer <- 6000
mask.spacing <- 50
trap.spacing <- 500

## True parameter values.
D <- 0.0416
g0 <- 1
sigma <- 1250
kappa <- 70
alpha <- 750 # Arbitrarily chosen value.
detectpars <- list(g0 = g0, sigma = sigma)

## Setting up mask and traps.
traps <- make.grid(nx = 3, ny = 1, spacing = trap.spacing, detector = "proximity")
ntraps <- nrow(traps)
mask <- make.mask(traps, spacing = mask.spacing, type = "trapbuffer", buffer = buffer)
nmask <- nrow(mask)
A <- attr(mask, "area")
res <- matrix(0, nrow = nsims, ncol = 5)
colnames(res) <- c("D", "sigma", "kappa", "alpha", "maxgrad")

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
  capthist.ang <- capthist
  capthist.ang[capthist.ang == 1] <- radians[capthist.ang == 1]
  captdists <- distances(as.matrix(detections), as.matrix(traps))
  distests <- function(x, alpha){
    betas <- alpha/x
    c(rgamma(1, shape = alpha, rate = betas[1]),
      rgamma(1, shape = alpha, rate = betas[2]),
      rgamma(1, shape = alpha, rate = betas[3]))
  }
  estdists <- t(apply(captdists, 1, distests, alpha = alpha))
  estdists <- array(estdists, dim = c(dim(estdists)[1], 1, dim(captdists)[2]))
  capthist.dist <- as.array(estdists*capthist)
  capthist.joint <- array(0, dim = c(dim(capthist), 2))
  capthist.joint[, , , 1] <- capthist.ang
  capthist.joint[, , , 2] <- capthist.dist
  fit.joint <- admbsecr(capthist.joint, traps = traps, mask = mask, method = "angdist",
                        detfn = "hn", bounds = list(alpha = c(0, 10000)),
                        sv = c(D = D, sigma = sigma, alpha = alpha, kappa = kappa),                     
                        fix = list(g0 = 1), scalefactors = c(D = 100000, kappa = 15))
  res[i, ] <- c(coef(fit.joint), fit.joint$maxgrad)
}
