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

p <- c(log(0.1125153), logit(0.95), log(750), log(10))
## Carrying out angle analysis with admbsecr(). Start values same as
## above provided, need to change to sv = sv. Or else use sv = "auto".
sv <- c("D" = exp(p[1]), "g0" = invlogit(p[2]), "sigma" = exp(p[3]), "kappa" = exp(p[4]))
anglefit <- admbsecr(capt = radians, traps = traps, mask = mask, #fix = list(g0 = 1),
                     sv = "auto", admbwd = admb.dir, method = "ang")
