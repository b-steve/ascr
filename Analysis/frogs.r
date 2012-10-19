## Letting R know where everything is.
admbsecr.dir <- "~/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Analysis", sep = sep)
func.dir <- paste(admbsecr.dir, "R", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", sep = sep)

## Running setup code.
setwd(work.dir)
source("frogsetup.r")

## Carrying out simple SECR analysis.

## With secr package.
simplefit1 <- secr.fit(capt, model = list(D~1, g0~1, sigma~1), mask = mask, verify = FALSE)
## With admbsecr().
simplefit2 <- admbsecr(capt, traps = traps, mask, sv = "auto", admbwd = admb.dir,
                       method = "simple")

## Carrying out TOA analysis.

## Create ssq matrix in advance.
ssqtoa <- apply(capt.toa, 1, toa.ssq, dists = dists)
## Starting value for TOA measurement standard deviation.
sigma.toa <- 0.0025
## Adding sigma.toa to estimates found above in the simple SECR for starting values.
start.beta <- c(simplefit1$fit$estimate, log(sigma.toa))
## Fitting with admbsecr(). Doesn't require start values, but can use sv = start.beta.
toafit1 <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto", 
                    admbwd = admb.dir, method = "toa")

## Carrying out signal strength analysis.

## Fitting with secr.fit().
ssfit1 <- secr.fit(sscapt, model = list(D~1, g0~1, sigma~1), detectfn = 10, mask = mask,
                   verify = FALSE, steptol = 1e-4)
## Fitting with admbsecr().
ssfit2 <- admbsecr(capt.ss, traps = traps, mask, sv = "auto", cutoff = 150,
                   admbwd = admb.dir, method = "ss")

## Carrying out analysis with both signal strength and TOA information incorporated.
## Only possible with admbsecr().
jointfit <- admbsecr(capt = capt.joint, traps = traps, mask = mask, sv = "auto",
                     cutoff = 150, admbwd = admb.dir, method = "sstoa")
