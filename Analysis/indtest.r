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
D <- 1750
g0 <- 1
sigma <- 8
pars <- c(D = D, g0 = g0, sigma = sigma)
bounds <- NULL
set.seed(seed)
## Calls per frog
cpf <- 10

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")

res <- list()
res.se <- list()
for (i in 1:nsims){
  capt <- sim.capt(traps = traps, calls = cpf, mask = mask, pars = pars)
  fit <- admbsecr(capt, traps = traps, mask = mask, sv = pars, fix = list(g0 = 1))
  fit.se <- se.correct(fit, calls = cpf, size = 100)
  res[[i]] <- fit
  res.se[[i]] <- fit.se
  save.image(file = "indtest.RData")
  print(c(date(), i))
}

save.image(file = "indtest.RData")

load(file = "~/admbsecr/Analysis/indtest.RData")

library(plyr)
library(admbsecr)

pars <- laply(res, coef)
pars.c <- laply(res.se, function(object) object$se.correct$coefficients.corrected)
ses <- laply(res, stdEr)
ses.c <- laply(res.se, function(object) object$se.correct$se.corrected)

Ds <- pars[, 1]
Ds.c <- pars.c[, 1]
D.se <- ses[, 1]
D.se.c <- ses.c[, 1]
