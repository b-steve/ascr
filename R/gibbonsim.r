library("secr")
library("CircStats")
library("inline")
library("Rcpp")

if (.Platform$OS == "unix"){
    source("/home/ben/SECR/R/helpers.r")
    source("/home/ben/SECR/R/admbsecr.r")
    load("/home/ben/SECR/Data/Gibbons/gibbons_data.RData")
    admb.dir <- "/home/ben/SECR/ADMB"
    dat.dir <- "/home/ben/SECR/Data/Gibbons/gibbons.txt"
} else if (.Platform$OS == "windows"){
    source("C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R\\helpers.r")
    source("C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R\\admbsecr.r")
    load("C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Gibbons\\gibbons_data.RData")
    admb.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\ADMB"
    dat.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Gibbons\\gibbons.txt"
}

nsims <- 4000
buffer <- 6000
mask.spacing <- 100
trap.spacing <- 2500

## True parameter values:
D <- 0.0416
g0 <- 0.99999
sigma <- 1250
kappa <- 140
truepars <- c(D = D, g0 = g0, sigma = sigma, kappa = kappa)
detectpars <- list(g0 = g0, sigma = sigma)

## Setting up mask and traps
traps <- make.grid(nx = 3, ny = 1, spacing = trap.spacing, detector = "proximity")
ntraps <- nrow(traps)
mask <- make.mask(traps, spacing = mask.spacing, type = "trapbuffer", buffer = buffer)
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances.cpp(as.matrix(traps), as.matrix(mask))
mask.angs <- angles.cpp(as.matrix(traps), as.matrix(mask))

simpleres <- matrix(0, nrow = nsims, ncol = 3)
angres <- matrix(0, nrow = nsims, ncol = 4)
colnames(simpleres) <- c("D", "g0", "sigma")
colnames(angres) <- c("D", "g0", "sigma", "kappa")
print(c("start", date()))
for (i in 1:nsims){
  if (i %% 100 == 0){
    print(c(i, date()))
  }
  ## Simulating data and setting things up for analysis
  popn <- sim.popn(D = D, core = traps, buffer = buffer)
  capthist <- sim.capthist(traps, popn, detectfn = 0, detectpar = detectpars, noccasions = 1,
                           renumber = FALSE)
  n <- nrow(capthist)
  ndets <- sum(capthist)
  cue.ids <- unique(as.numeric(rownames(capthist)))
  detections <- popn[cue.ids, ]
  radians <- t(angles.cpp(as.matrix(traps), as.matrix(detections)))
  radians <- array(radians, c(dim(radians), 1))
  errors <- array(rvm(ntraps*n, mean = 0, k = kappa), dim(radians))
  radians <- (radians + errors) %% (2*pi)
  radians[radians == 0] <- 2*pi
  radhist <- capthist
  radhist[radhist == 1] <- radians[radhist == 1]
  
  p <- c(log(D), invlogit(g0), log(sigma))
  ## Straightforward secr model using admbsecr()
  simplefit <- admbsecr(capt = capthist, traps = traps, mask = mask, sv = truepars[1:3],
                        admbwd = admb.dir, method = "simple", verbose = FALSE, autogen = FALSE)
  angfit <- admbsecr(capt = radhist, traps = traps, mask = mask, sv = truepars,
                     angs = mask.angs, admbwd = admb.dir, method = "ang", verbose = FALSE,
                     autogen = FALSE)
  simpleres[i, ] <- coef(simplefit)
  angres[i, ] <- coef(angfit)
}
print(c("end", date()))


plot(simpleres[, 1], rep(1, 10), ylim = c(0, 3))
points(angres[, 1], rep(2, 10))
abline(v = truepars[1], lty = "dotted")
for (i in 1:10){
  lines(x = c(simpleres[i, 1], angres[i, 1]), y = 1:2)
}
