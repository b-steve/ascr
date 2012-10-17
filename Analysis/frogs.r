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
## Fitting with nlm().
toafit1 <- nlm(f = secrlikelihood.toa1, p = start.beta, capthist=capt.toa,
                      mask = mask, dists = dists, ssqtoa = ssqtoa, trace = TRUE)
## Fitting with admbsecr(). Doesn't require start values.
toafit2 <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto", 
                    admbwd = admb.dir, method = "toa")

## Carrying out signal strength analysis.

## Fitting with nlm().
startval <- c(log(7000), 190, 0, log(6))
ssfit1 <- nlm(f = secrlikelihood.ss, p = startval, capthist = capt.ss, mask = mask,
              dists = dists, cutoff = 150, trace = TRUE)
ests <- c(exp(ssfit1$estimate[1]), ssfit1$estimate[2:3], exp(ssfit1$estimate[4]))
## Fitting with secr.fit().
ssfit2 <- secr.fit(sscapt, model = list(D~1, g0~1, sigma~1), detectfn = 10, mask = mask,
                   verify = FALSE, steptol = 1e-4)
## Fitting with admbsecr().
ssfit3 <- admbsecr(capt.ss, traps = traps, mask, sv = "auto", cutoff = 150,
                   admbwd = admb.dir, method = "ss")

## Carrying out analysis with both signal strength and TOA information incorporated.
## Only possible with admbsecr().
jointfit <- admbsecr(capt = capt.joint, traps = traps, mask = mask, sv = "auto",
                     cutoff = 150, admbwd = admb.dir, method = "sstoa")
