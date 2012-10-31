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
## Setup for simulations.
nsims <- 1
buffer <- 35
mask.spacing <- 50

## True parameter values.
seed <- 5253
D <- 4450
g0 <- 0.99999
sigma <- 5.60
sigmatoa <- 0.002
ssb0 <- 170
ssb1 <- -2.50
sigmass <- 7.00
cutoff <- 150
truepars <- c(D = D, g0 = g0, sigma = sigma, sigmatoa = sigmatoa,
              ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)

set.seed(seed)
## Inverse of speed of sound (in ms per metre).
invsspd <- 1000/330

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(traps), as.matrix(mask))
simprobs <- NULL
toaprobs <- NULL
ssprobs <- NULL
jointprobs <- NULL
simpleres <- matrix(0, nrow = nsims, ncol = 3)
toares <- matrix(0, nrow = nsims, ncol = 4)
ssres <- matrix(0, nrow = nsims, ncol = 4)
jointres <- matrix(0, nrow = nsims, ncol = 5)
colnames(simpleres) <- c("D", "g0", "sigma")
colnames(toares) <- c("D", "g0", "sigma", "sigmatoa")
colnames(ssres) <- c("D", "ssb0", "ssb1", "sigmass")
colnames(jointres) <- c("D", "sigmatoa", "ssb0", "ssb1", "sigmass")

## Carrying out simulation.
for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else if (i %% 100 == 0){
    print(c(i, date()))
  }
  ## Simulating data and setting things up for analysis.
  popn <- sim.popn(D = D, core = traps, buffer = buffer)
  capthist <- sim.capthist(traps, popn, detectfn = 10, detectpar = detectpars,
                       noccasions = 1, renumber = FALSE)
  n <- nrow(capthist)
  ndets <- sum(capthist)
  ## IDs for detected animals.
  cue.ids <- unique(as.numeric(rownames(capthist)))
  ## Cartesian coordinates of detected animals.
  detections <- popn[cue.ids, ]
  ## Distances from detected animals to traps.
  dists <- t(distances(as.matrix(traps), as.matrix(detections)))
  capthist.ss <- array(0, dim = dim(capthist))
  capthist.ss[capthist == 1] <- attr(capthist, "signalframe")[, 1]
  ## Generating times of arrival.
  ## Time of call itself doesn't provide any information, this comes
  ## solely from the *differences* in arrival times between traps. We
  ## can assume that animal with ID = t calls at time t without loss
  ## of generality.
  capthist.toa <- array(0, dim = dim(capthist))
  for (j in 1:n){
    for (k in 1:ntraps){
      if (capthist[j, 1, k] == 1){
        dist <- dists[j, k]
        meantoa <- cue.ids[j] + invsspd*dist/1000
        capthist.toa[j, 1, k] <- rnorm(1, meantoa, sigmatoa)
      } else {
        capthist.toa[j, 1, k] <- 0
      }
    }
  }
  capthist.joint <- array(c(capthist.ss, capthist.toa), dim = c(dim(capthist), 2))
  ## Straightforward SECR model using admbsecr().
  simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                            sv = truepars[1:3], admbwd = admb.dir,
                            method = "simple", verbose = FALSE), silent = TRUE)
  if (class(simplefit) == "try-error"){
    simplefit <- try(admbsecr(capt = capthist, traps = traps, mask = mask,
                              sv = "auto", admbwd = admb.dir, method = "simple",
                              verbose = FALSE), silent = TRUE)
  }
  if (class(simplefit) == "try-error"){
    simplecoef <- NA
    simprobs <- c(simprobs, i)
  } else {
    simplecoef <- coef(simplefit)
  }
  ## SECR model using supplementary TOA data.
  ssqtoa <- apply(capthist.toa, 1, toa.ssq, dists = mask.dists)
  toafit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask,
                         sv = truepars[1:4], ssqtoa = ssqtoa, admbwd = admb.dir,
                         method = "toa", verbose = FALSE), silent = TRUE)
  if (class(toafit) == "try-error"){
    toafit <- try(admbsecr(capt = capthist.toa, traps = traps, mask = mask,
                           sv = "auto", ssqtoa = ssqtoa, admbwd = admb.dir,
                           method = "toa", verbose = FALSE), silent = TRUE)
  }
  if (class(toafit) == "try-error"){
    toacoef <- NA
    toaprobs <- c(toaprobs, i)
  } else {
    toacoef <- coef(toafit)
  }
  ## SECR model using supplementary signal strength data.
  ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                        sv = truepars[c(1, 5:7)], cutoff = cutoff,
                        admbwd = admb.dir, method = "ss", verbose = FALSE),
               silent = FALSE)
  if (class(ssfit) == "try-error"){
    ssfit <- try(admbsecr(capt = capthist.ss, traps = traps, mask = mask,
                          sv = "auto", cutoff = cutoff, admbwd = admb.dir,
                          method = "ss", verbose = FALSE), silent = FALSE)
  }
  if (class(ssfit) == "try-error"){
    sscoef <- NA
    ssprobs <- c(ssprobs, i)
  } else {
    sscoef <- coef(ssfit)
  }
  ## SECR model using joint TOA and signal strength data.
  jointfit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                           sv = truepars[c(1, 4:7)], cutoff = cutoff,
                           ssqtoa = ssqtoa, admbwd = admb.dir, method = "sstoa",
                           verbose = FALSE), silent = TRUE)
  if (class(jointfit) == "try-error"){
    jointfit <- try(admbsecr(capt = capthist.joint, traps = traps, mask = mask,
                             sv = "auto", cutoff = cutoff, ssqtoa = ssqtoa,
                             admbwd = admb.dir, method = "sstoa",
                             verbose = FALSE), silent = TRUE)
  }
  if (class(jointfit) == "try-error"){
    jointcoef <- NA
    jointprobs <- c(jointprobs, i)
  } else {
    jointcoef <- coef(jointfit)
  }
  simpleres[i, ] <- simplecoef
  toares[i, ] <- toacoef
  ssres[i, ] <- sscoef
  jointres[i, ] <- jointcoef
  if (i == nsims){
    print(c("end", date()))
  }
}

## To write the simulation results to a file.
## write.table(simpleres, "~/admbsecr/Results/frogs/1/simpleres.txt", row.names = FALSE)
## write.table(toares, "~/admbsecr/Results/frogs/1/toares.txt", row.names = FALSE)
## write.table(ssres, "~/admbsecr/Results/frogs/1/ssres.txt", row.names = FALSE)
## write.table(jointres, "~/admbsecr/Results/frogs/1/jointres.txt", row.names = FALSE)

## To read in simulation results from a file.
resfile <- "/home/ben/admbsecr/Results/frogs/1/"
source(paste(resfile, "pars.r", sep = ""))
simpleres <- read.table(paste(resfile, "simpleres.txt", sep = ""), header = TRUE)
toares <- read.table(paste(resfile, "toares.txt", sep = ""), header = TRUE)
ssres <- read.table(paste(resfile, "ssres.txt", sep = ""), header = TRUE)
jointres <- read.table(paste(resfile, "jointres.txt", sep = ""), header = TRUE)

## Assigning the columns to vectors.
for (i in colnames(simpleres)){
  name <- paste("sim", i, sep = "")
  assign(name, simpleres[, i])
}
for (i in colnames(toares)){
  name <- paste("toa", i, sep = "")
  assign(name, toares[, i])
}
for (i in colnames(ssres)){
  name <- paste("ss", i, sep = "")
  assign(name, ssres[, i])
}
for (i in colnames(jointres)){
  name <- paste("joint", i, sep = "")
  assign(name, jointres[, i])
}



## Two different bandwidth selections.
dsimD <- density(simD)
dtoaD <- density(toaD)
dssD <- density(ssD)
djointD <- density(jointD)
xs <- c(dsimD$x, dtoaD$x, dssD$x, djointD$x)
ys <- c(dsimD$y, dtoaD$y, dssD$y, djointD$y)



##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(dsimD, col = "blue")
lines(dtoaD, col = "red")
lines(dssD, col = "green")
lines(djointD, col = "purple")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density",
      xlab = expression(hat(D)), ylab = "Density")
legend(x = "topright", legend = c("SECR", "SECR + TOA", "SECR + SS", "SECR + TOA + SS"),
       col = c("blue", "red", "green", "purple"), lty = 1)
##dev.off()
