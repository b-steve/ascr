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
datasource <- "Res"
source("frogsetup.r")
cutoff <- ifelse(datasource == "Res", 130, 150)

cpi <- c(rep(17, 4), rep(16, 2), rep(15, 2))

## Carrying out simple SECR analysis.

## With secr package.
##simplefit1 <- secr.fit(capt, model = list(D~1, g0~1, sigma~1), mask = mask, verify = FALSE, method = "BFGS")

## With admbsecr().
simplefit.hn <- admbsecr(capt, traps = traps, mask, sv = "auto", method = "simple")
simplefit.th <- admbsecr(capt, traps = traps, mask, sv = c(shape = 2.63, scale = 3),
                         method = "simple", detfn = "th")
##simplefit.logth <- admbsecr(capt, traps = traps, mask,
##                            sv = c(shape1 = 17.1493, shape2 = 2.9919, scale = -0.0204),
##                            method = "simple", detfn = "logth")
simplefit.hr1 <- admbsecr(capt, traps = traps, mask, sv = c(sigma = 7, z = 12),
                          method = "simple", detfn = "hr")
simplefit.hr2 <- admbsecr(capt, traps = traps, mask, sv = c(sigma = 7, z = 12),
                          fix = list(g0 = 1), method = "simple", detfn = "hr")

## Carrying out TOA analysis.

## Fitting with admbsecr(). Doesn't require start values.
toafit.hn <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto",
                      method = "toa")
toafit.th <- admbsecr(capt = capt.toa, traps = traps, mask = mask,
                      sv = c(shape = 2.63, scale = 3),
                      method = "toa", detfn = "th")
##toafit.logth <- admbsecr(capt = capt.toa, traps = traps, mask = mask,
##                         sv = c(shape1 = 17.1493, shape2 = 2.9919, scale = -0.0204),
##                         bounds = list(shape1 = c(0, 10), scale = c(-5, -0.01)),
##                         method = "toa", detfn = "logth")
if (datasource == "Res") hr1.bounds <- list(D = c(0, 10000)) else hr1.bounds <- NULL
toafit.hr1 <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto",
                       bounds = hr1.bounds, method = "toa", detfn = "hr")
if (datasource == "Res") hr2.bounds <- list(sigma = c(0, 1000)) else hr2.bounds <- NULL
toafit.hr2 <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto",
                       bounds = hr2.bounds, fix = list(g0 = 1), method = "toa", detfn = "hr")

## Carrying out signal strength analysis.

## Fitting with secr.fit().
##ssfit1 <- secr.fit(sscapt, model = list(D~1, g0~1, sigma~1), detectfn = 10, mask = mask,
##                   verify = FALSE, steptol = 1e-4)

## Fitting with admbsecr().
ssfit2 <- admbsecr(capt.ss, traps = traps, mask, cutoff = cutoff, sv = "auto",
                   admbwd = admb.dir, method = "ss")
ssfit2.log <- admbsecr(capt.ss, traps = traps, mask, cutoff = cutoff, sv = "auto",
                       bounds = list(ssb1 = c(0, 5)),
                       admbwd = admb.dir, method = "ss", detfn = "log")

## Carrying out analysis with both signal strength and TOA information incorporated.
## Only possible with admbsecr().
jointfit <- admbsecr(capt = capt.joint, traps = traps, mask = mask,
                     bounds = list(ssb0 = c(cutoff, 1e8)),
                     cutoff = cutoff, method = "sstoa")
jointfit.log <- admbsecr(capt = capt.joint, traps = traps, mask = mask,
                         cutoff = cutoff, method = "sstoa", detfn = "log")


## Profile likelihood confidence interval? Doesn't work properly.
jointfit.prof <- admbsecr(capt = capt.joint, traps = traps, mask = mask,
                     bounds = list(ssb0 = c(cutoff, 1e8)),
                     cutoff = cutoff, admbwd = admb.dir, method = "sstoa", profpars = "D")


## Analysis from SS/TOA frog paper (with cpi):
simplefit.hn <- admbsecr(capt, traps = traps, mask, sv = "auto", method = "simple",
                         cpi = cpi)
simplefit.th <- admbsecr(capt, traps = traps, mask, sv = c(shape = 2.63, scale = 3),
                         method = "simple", detfn = "th", cpi = cpi)
simplefit.hr <- admbsecr(capt, traps = traps, mask, sv = c(sigma = 7, z = 12),
                         fix = list(g0 = 1), method = "simple", detfn = "hr",
                         cpi = cpi)
simplefit.hn.c <- se.correct(simplefit.hn, 100)
simplefit.th.c <- se.correct(simplefit.th, 100)
simplefit.hr.c <- se.correct(simplefit.hr, 100)

## Fitting with admbsecr(). Doesn't require start values.
toafit.hn <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto",
                      method = "toa", cpi = cpi)
toafit.th <- admbsecr(capt = capt.toa, traps = traps, mask = mask,
                      sv = c(shape = 2.63, scale = 3),
                      method = "toa", detfn = "th", cpi = cpi)
if (datasource == "Res") hr2.bounds <- list(sigma = c(0, 1000)) else hr2.bounds <- NULL
toafit.hr2 <- admbsecr(capt = capt.toa, traps = traps, mask = mask, sv = "auto",
                       bounds = hr2.bounds, fix = list(g0 = 1), method = "toa", detfn = "hr",
                       cpi = cpi)
toafit.hn.c <- secorrect(toafit.hn, 100)
toafit.th.c <- secorrect(toafit.th, 100)
toafit.hr.c <- secorrect(toafit.hr, 100)

ssfit <- admbsecr(capt.ss, traps = traps, mask, cutoff = cutoff, sv = "auto",
                   method = "ss", cpi = cpi)
ssfit.log <- admbsecr(capt.ss, traps = traps, mask, cutoff = cutoff, sv = "auto",
                       bounds = list(ssb1 = c(0, 5)),  method = "ss", detfn = "log",
                       cpi = cpi)
ssfit.c <- secorrect(ssfit, 100)
ssfit.log.c <- secorrect(ssfit.log)

jointfit <- admbsecr(capt = capt.joint, traps = traps, mask = mask,
                     bounds = list(ssb0 = c(cutoff, 1e8)),
                     cutoff = cutoff, admbwd = admb.dir, method = "sstoa",
                     cpi = cpi)
jointfit.log <- admbsecr(capt = capt.joint, traps = traps, mask = mask,
                         cutoff = cutoff, method = "sstoa", detfn = "log",
                         cpi = cpi)
jointfit.c <- secorrect(ssfit, 100)
jointfit.log.c <- secorrect(jointfit.log, 100)
