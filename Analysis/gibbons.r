## Letting R know where everything is.
admbsecr.dir <- "~/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
## admb.dir only required when user provides the .tpl file.
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Analysis", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", sep = sep)

## Load admbsecr either using devtools or as a library.
setwd(admbsecr.dir)
library(devtools)
load_all()
##library(admbsecr)

## Running setup code.
setwd(work.dir)
library(secr)
source("gibbonsetup.r")

simplefit.hn <- admbsecr(capt = capthist, traps = traps, mask = mask, fix = list(g0 = 1),
                         sv = "auto", admbwd = admb.dir, method = "simple")
simplefit.th <- admbsecr(capt = capthist, traps = traps, mask = mask,
                         sv = c(shape = -5, scale = -0.005), bounds = list(shape = c(-15, 0)),
                         detfn = "th", method = "simple")
simplefit.hr1 <- admbsecr(capt = capthist, traps = traps, mask = mask, sv = "auto",
                          detfn = "hr", method = "simple", clean = FALSE)
simplefit.hr2 <- admbsecr(capt = capthist, traps = traps, mask = mask, fix = list(g0 = 1),
                          sv = c(sigma = 1000, z = 5.5), detfn = "hr", method = "simple",
                          bounds = list(sigma = c(0, 4000)), clean = FALSE, trace = TRUE)

## Carrying out angle analysis with admbsecr().
anglefit.hn <- admbsecr(capt = radians, traps = traps, mask = mask, fix = list(g0 = 1),
                        sv = "auto", method = "ang")
anglefit.th <- admbsecr(capt = radians, traps = traps, mask = mask,
                        sv = c(shape = -5, scale = -0.005), bounds = list(shape = c(-15, 0)),
                        detfn = "th", method = "ang")
anglefit.hr1 <- admbsecr(capt = radians, traps = traps, mask = mask, sv = "auto",
                         detfn = "hr", method = "ang")
anglefit.hr2 <- admbsecr(capt = radians, traps = traps, mask = mask, fix = list(g0 = 1),
                         sv = "auto", detfn = "hr", method = "ang")

## Comparison of fitted detection functions.
hn <- function(x, coef){
  g0 <- coef[2]
  sigma <- coef[3]
  g0*exp(-x^2/(2*sigma^2))
}
ss <- function(x, c, coef){
  ssb0 <- coef[2]
  ssb1 <- coef[3]
  sigmass <- coef[4]
  1 - pnorm(c, ssb0 + ssb1*x, sigmass)
}
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
th <- function(x, coef){
  shape <- coef[2]
  scale <- coef[3]
  0.5 - 0.5*erf(shape - scale*x)
}
hr <- function(x, coef){
  g0 <- coef[2]
  sigma <- coef[3]
  z <- coef[4]
  g0*(1 - exp(-(x/sigma)^(-z)))
}

## Comparison using simple method.
x <- seq(0, 3000, 1)
plot(x, hn(x, c(0, 1, coef(simplefit.hn)[2])), type = "l", ylim = 0:1)
lines(x, th(x, coef(simplefit.th)), col = "green")
lines(x, hr(x, coef(simplefit.hr1)), col = "blue")
hr2pars <- c(0, 1, coef(simplefit.hr2)[2:3])
lines(x, hr(x, hr2pars), col = "brown")

## Comparison using angles method.
x <- seq(0, 4000, 1)
plot(x, hn(x, c(0, 1, coef(anglefit.hn)[2])), type = "l", ylim = 0:1)
lines(x, th(x, coef(anglefit.th)), col = "green")
lines(x, hr(x, coef(anglefit.hr1)), col = "blue")
hr2pars <- c(0, 1, coef(anglefit.hr2)[2:3])
lines(x, hr(x, hr2pars), col = "brown")
