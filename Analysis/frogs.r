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
source("frogsetup.r")

## Carrying out simple SECR analysis.

## With secr package.
simplefit1 <- secr.fit(capt, model = list(D~1, g0~1, sigma~1), mask = mask, verify = FALSE)
## With admbsecr().
simplefit2 <- admbsecr(capt, traps = traps, mask, sv = "auto", admbwd = admb.dir,
                       method = "simple")

## Carrying out TOA analysis.

## Fitting with admbsecr(). Doesn't require start values.
toafit1 <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto",
                    admbwd = admb.dir, method = "toa")

## Carrying out signal strength analysis.

## Fitting with secr.fit().
ssfit1 <- secr.fit(sscapt, model = list(D~1, g0~1, sigma~1), detectfn = 10, mask = mask,
                   verify = FALSE, steptol = 1e-4)
## Fitting with admbsecr().
ssfit2 <- admbsecr(capt.ss, traps = traps, mask, cutoff = 150, sv = "auto",
                   admbwd = admb.dir, method = "ss", autogen = FALSE)

## Carrying out analysis with both signal strength and TOA information incorporated.
## Only possible with admbsecr().
##sv <- c(coef(toafit1)[c(1, 4)], coef(ssfit2)[2:4])
jointfit <- admbsecr(capt = capt.joint, traps = traps, mask = mask,
                     bounds = list(D = c(0, 1e5), ssb1 = c(-1, 0)),
                     cutoff = 150, admbwd = admb.dir, method = "sstoa",
                     autogen = FALSE, trace = TRUE)
