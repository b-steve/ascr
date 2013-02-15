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

## Load admbsecr either using devtools or as a library.
setwd(admbsecr.dir)
library(devtools)
load_all()
##library(admbsecr)

## Running setup code.
setwd(work.dir)
library(secr)
source("whalesetup.r")

distfit87 <- admbsecr(capthist87.dist, traps = real.traps, mask = mask87,
                      sv = c(D = 10, g0 = 0.95, sigma = 100, alpha = 5),
                      method = "dist")

distfit01.hn <- admbsecr(capthist01.dist, traps = real.traps, mask = mask01,
                      sv = c(D = 10, g0 = 0.95, sigma = 100, alpha = 2),
                      method = "dist")
distfit01.th <- admbsecr(capthist01.dist, traps = real.traps, mask = mask01,
                         sv = c(D = 10, shape = -1, scale = -0.01, alpha = 2),
                         method = "dist", detfn = "th")
distfit01.hr <- admbsecr(capthist01.dist, traps = real.traps, mask = mask01,
                         sv = c(D = 10, g0 = 0.95, sigma = 100, z = 5, alpha = 2),
                         method = "dist", detfn = "hr")

mrdsfit01.hn <- admbsecr(capthist.mrds, traps = real.traps, mask = mask01,
                      sv = c(D = 10, g0 = 0.95, sigma = 100),
                      method = "mrds")
mrdsfit01.th <- admbsecr(capthist.mrds, traps = real.traps, mask = mask01,
                         sv = c(D = 5, shape = -0.24, scale = -0.003),
                         bounds = list(D = c(0, 10)),
                         method = "mrds", detfn = "th", trace = TRUE)
mrdsfit01.hr <- admbsecr(capthist.mrds, traps = real.traps, mask = mask01,
                         sv = c(D = 10, g0 = 0.95, sigma = 100, z = 5),
                         method = "mrds", detfn = "hr")

## dist fit with separate detection functions for each trap:
distfit01tc <- disttrapcov(capt = capthist01.dist, mask = mask01, traps = real.traps,
                           sv = c(10, 0.5, 100, 0.5, 100, 5), admb.dir = admb.dir,
                           clean = TRUE, verbose = FALSE, trace = FALSE)

## mrds fit with seperate detection functions for each trap:
mrdsfit01tc <- mrdstrapcov(capt = capthist.mrds, mask = mask01, traps = real.traps,
                           sv = c(10, 0.5, 100, 0.5, 100), admb.dir = admb.dir,
                           clean = TRUE, verbose = FALSE, trace = FALSE)

## Comparing the models to see if the extra g0 and sigma parameters are required:

## Likelihood-ratio tests both strongly suggest g01 != g02, sigma1 != sigma2
LRTS <- 2*(logLik(mrdsfit01tc) - logLik(mrdsfit01))
1 - pchisq(LRTS, 2)
LRTS <- 2*(logLik(distfit01tc) - logLik(distfit01))
1 - pchisq(LRTS, 2)

## AICs also show the same thing.
AIC(mrdsfit01tc)
AIC(mrdsfit01)
AIC(distfit01tc)
AIC(distfit01)
## Note that you cannot compare mrds and dist fits using AIC here as
## the same data is not used (distances for whales seen by both
## observers are averaged for mrds).

## Contour plot example:
contours(distfit01tc, 1, ylim = c(-10, 1000), xlim = c(-1000, 1000))

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
x <- seq(0, 1000, 1)
## Comparison of detection functions without TOA.
plot(x, hn(x, coef(distfit01.hn)), type = "l", ylim = 0:1)
lines(x, th(x, coef(distfit01.th)), col = "green")
lines(x, hr(x, coef(distfit01.hr)), col = "blue")

plot(x, hn(x, coef(mrdsfit01.hn)), type = "l", ylim = 0:1)
lines(x, th(x, coef(mrdsfit01.th)), col = "green")
lines(x, hr(x, coef(mrdsfit01.hr)), col = "blue")
