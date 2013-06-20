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
work.dir <- paste(admbsecr.dir, "Analysis/tests", sep = sep)
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
D <- 1750
shape <- 2.39
scale <- 5.35
pars <- c(D = D, shape = shape, scale = scale)
bounds <- NULL
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

res <- list()
res.se <- list()
for (i in 1:nsims){
  capt <- sim.capt(traps = traps, calls = cpf, mask = mask, pars = pars, detfn = "th")
  system.time({
  fit <- try.admbsecr(sv = pars, capt = capt, traps = traps, mask = mask,
                      detfn = "th", bounds = bounds, scalefactors = scalefactors)})
  if (class(fit)[1] == "try-error"){
    res[[i]] <- "error"
  } else {
    fit.se <- se.correct(fit, calls = cpf, size = 100)
    res[[i]] <- fit
    res.se[[i]] <- fit.se
  }
  save.image(file = "indtest.RData")
  print(c(date(), i))
}

save.image(file = "indtest.RData")

load(file = "~/admbsecr/Results/ind/1")

library(plyr)
library(admbsecr)

pars <- laply(res, coef)
pars.c <- laply(res.se, function(object) object$se.correct$coefficients.corrected)
ses <- laply(res, stdEr)
ses.c <- laply(res.se, function(object) object$se.correct$se.corrected)

ns <- laply(res.se, function(object) dim(object$se.correct$boots)[1])

Ds <- pars[, 1]
Ds.c <- pars.c[, 1]
D.se <- ses[, 1]
D.se.c <- ses.c[, 1]

cis <- cis.c <- matrix(0, nrow = nsims, ncol = 2)
for (i in 1:nsims){
  cis[i, ] <- Ds[i] + c(-1, 1)*1.96*D.se[i]
  cis.c[i, ] <- Ds.c[i] + c(-1, 1)*1.96*D.se.c[i]
}
cov <- mean(apply(cis, 1, function(x, mu) mu >= x[1] & mu <= x[2], mu = D))
cov.c <- mean(apply(cis.c, 1, function(x, mu) mu >= x[1] & mu <= x[2], mu = D))
