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

distfit01 <- admbsecr(capthist01.dist, traps = real.traps, mask = mask01,
                      sv = c(D = 10, g0 = 0.95, sigma = 100, alpha = 2),
                      method = "dist")

mrdsfit01 <- admbsecr(capthist.mrds, traps = real.traps, mask = mask01,
                      sv = c(D = 10, g0 = 0.95, sigma = 100),
                      method = "mrds")
