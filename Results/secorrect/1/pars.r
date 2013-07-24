wd <- getwd()

## Setup for simulations.
seed <- 15338
set.seed(seed)
nsims <- 4
nboots <- 100
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
setwd("~/admbsecr/Data")
mics <- read.csv(file = "array1a-top.csv")
micnames <- 1:dim(mics)[1]
mics <- cbind(micnames, mics)
traps <- read.traps(data = mics, detector = "signal")
setwd(wd)
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
