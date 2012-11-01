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
library(plyr)

## Loading the admbsecr library.
setwd(admbsecr.dir)
library(devtools)
load_all(".")
##library(admbsecr)

## Creating trap objects.
fake.traps <- data.frame(name = 1:2, x = c(1e-3, -1e-3), y = c(-1e-3, 1e-3))
real.traps <- data.frame(name = 1:2, x = rep(0, 2), y = rep(0, 2))
fake.traps <- read.traps(data = fake.traps, detector = "proximity")
options(warn = -1)
real.traps <- read.traps(data = real.traps, detector = "proximity")
options(warn = 1)

setwd(work.dir)
## setup for simulations.
nsims <- 500
buffer <- 3000
mask.spacing <- 45

## True parameter values.
seed <- 6843
D <- 7.10
g0 <- 0.29
sigma <- 245
alpha <- 5.6
truepars <- c(D = D, g0 = g0, sigma = sigma, alpha = alpha)
detectpars <- list(g0 = g0, sigma = sigma)

set.seed(seed)

## Setting up mask and traps.
ntraps <- nrow(traps)
mask <- make.mask(fake.traps, buffer = buffer, spacing = mask.spacing,
                  type = "trapbuffer")
nmask <- nrow(mask)
A <- attr(mask, "area")
mask.dists <- distances(as.matrix(real.traps), as.matrix(mask))
probs <- NULL
res <- matrix(0, nrow = nsims, ncol = 4)
colnames(res) <- c("D", "g0", "sigma", "alpha")

## Function for generating random gamma values:
estdists <- function(x, alpha){
  betas <- alpha/x
  c(rgamma(1, shape = alpha, rate = betas[1]),
    rgamma(1, shape = alpha, rate = betas[2]))
}

## Carrying out simulation.
for (i in 1:nsims){
  if (i == 1){
    print(c("start", date()))
  } else {
    print(c(i, date()))
  }
  ## Simulating data and setting things up for analysis.
  popn <- sim.popn(D = D, core = real.traps, buffer = buffer)

  capthist <- sim.capthist(real.traps, popn, detectfn = 0,
                           detectpar = detectpars, noccasions = 1,
                           renumber = FALSE)
  n <- nrow(capthist)
  ndets <- sum(capthist)
  cue.ids <- unique(as.numeric(rownames(capthist)))
  detections <- popn[cue.ids, ]
  dists <- distances(as.matrix(detections), as.matrix(real.traps))
  dists <- aaply(dists, 1, estdists, alpha = alpha)
  dists <- array(dists, dim = c(dim(dists)[1], 1, dim(dists)[2]))
  capthist.dist <- as.array(dists*capthist)
  fit <- try(admbsecr(capthist.dist, traps = real.traps, mask = mask,
                  sv = c(D = D, g0 = g0, sigma = sigma, alpha = alpha),
                  admbwd = admb.dir, method = "dist", autogen = FALSE),
             silent = TRUE)
  if (class(fit) == "try-error"){
    coef <- NA
    probs <- c(probs, i)
  }
  coef <- coef(fit)
  res[i, ] <- coef
}

## To write the simulation results to a file.
##respath <- paste(admbsecr.dir, "Results/whales/1/res.txt", sep = "/")
##write.table(res, respath, row.names = FALSE)

## To read in simulation results from a file.
resfile <- "~/admbsecr/Results/whales/1/"
source(paste(resfile, "pars.r", sep = ""))
res <- read.table(paste(resfile, "res.txt", sep = ""), header = TRUE)

## Assigning the columns to vectors.
for (i in colnames(res)){
  name <- paste("res", i, sep = "")
  assign(name, res[, i])
}

## Calculating densities.
dD <- density(resD)
dg0 <- density(resg0)
dsigma <- density(ressigma)
dalpha <- density(resalpha)
xs <- dD$x
ys <- dD$y

##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(dD, col = "blue")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density",
      xlab = expression(hat(D)), ylab = "Density")
##dev.off()
