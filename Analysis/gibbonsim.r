## Letting R know where everything is.
admbsecr.dir <- "/home/ben/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Analysis", sep = sep)
func.dir <- paste(admbsecr.dir, "R", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", "Gibbons", sep = sep)

library("secr")
library("CircStats")
library("inline")
library("Rcpp")

## Set working directory to that with the functions.
setwd(func.dir)
## Get SECR functions.
source("admbsecr.r")
source("autofuns.r")
source("helpers.r")
source("lhoodfuns.r")
source("tplmake.r")

## Setup for simulations.
nsims <- 2
buffer <- 6000
mask.spacing <- 50
trap.spacing <- 500

## True parameter values.
D <- 0.0416
g0 <- 0.99999
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
mask.dists <- distances.cpp(as.matrix(traps), as.matrix(mask))
mask.angs <- angles.cpp(as.matrix(traps), as.matrix(mask))
simprobs <- NULL
angprobs <- NULL
simpleres <- matrix(0, nrow = nsims, ncol = 3)
angres <- matrix(0, nrow = nsims, ncol = 5)
aangres <- matrix(0, nrow = nsims, ncol = 5)
colnames(simpleres) <- c("D", "g0", "sigma")
colnames(angres) <- colnames(aangres) <- c("D", "g0", "sigma", "kappa", "logLik")
for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else if (i %% 100 == 0){
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
  ## Straightforward SECR model using admbsecr()
  simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, sv = truepars[1:3],
                        admbwd = admb.dir, method = "simple", verbose = FALSE, autogen = FALSE),
                   silent = TRUE)
  if (class(simplefit) == "try-error"){
    simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask, sv = "auto",
                              admbwd = admb.dir, method = "simple", verbose = FALSE,
                              autogen = FALSE), silent = TRUE)
  }
  if (class(simplefit) == "try-error"){
    simplecoef <- NA
    simprobs <- c(simprobs, i)
  } else {
    simplecoef <- coef(simplefit)
  }
  ## SECR model using supplementary angle data
  angfit <- try(admbsecr(capt = radhist, traps = traps, mask = mask, sv = truepars,
                     angs = mask.angs, admbwd = admb.dir, method = "ang", verbose = FALSE,
                     autogen = FALSE), silent = TRUE)
  if (class(angfit) == "try-error"){
    angfit <- try(admbsecr(capt = radhist, traps = traps, mask = mask, sv = "auto",
                       angs = mask.angs, admbwd = admb.dir, method = "ang", verbose = FALSE,
                       autogen = FALSE), silent = TRUE)
  }
  if (class(angfit) == "try-error"){
    angcoef <- NA
    angprobs <- c(angprobs, i)
  } else {
    angcoef <- c(coef(angfit), logLik(angfit))
  }
  simpleres[i, ] <- simplecoef
  angres[i, ] <- angcoef
  if (i == nsims){
    print(c("end", date()))
  }
}

## To write the simulation results to a file.
##write.table(angres, "/home/ben/SECR/Results/3/angres.txt", row.names = FALSE)
##write.table(simpleres, "/home/ben/SECR/Results/3/simpleres.txt", row.names = FALSE)

## To read in simulation results from a file.
##resfile <- "/home/ben/SECR/Results/3/"
##source(paste(resfile, "pars.r", sep = ""))
##angres <- read.table(paste(resfile, "angres.txt", sep = ""), header = TRUE)
##simpleres <- read.table(paste(resfile, "simpleres.txt", sep = ""), header = TRUE)

## Assigning the columns to vectors.
for (i in colnames(simpleres)){
  name <- paste("sim", i, sep = "")
  assign(name, simpleres[, i])
}
for (i in colnames(angres)){
  name <- paste("ang", i, sep = "")
  assign(name, angres[, i])
}

## Two different bandwidth selections.
##dsimD <- density(simD)
##dangD <- density(angD)
dsimD <- density(simD, bw = "bcv")
dangD <- density(angD, bw = "bcv")
xs <- c(dsimD$x, dangD$x)
ys <- c(dsimD$y, dangD$y)


##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(dsimD, col = "blue")
lines(dangD, col = "red")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density", xlab = expression(hat(D)), ylab = "Density")
legend(x = "topright", legend = c("SECR", "SECR + angles"), col = c("blue", "red"), lty = 1)
##dev.off()
