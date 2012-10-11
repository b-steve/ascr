library("secr")
library("Rcpp")
library("inline")

## Set working directory to that with the functions.
setwd(func.dir)
## Get SECR functions.
source("admbsecr.r")
source("autofuns.r")
source("helpers.r")
source("lhoodfuns.r")
source("tplmake.r")

setwd(dat.dir)
## Get and manipulate data.
mics <- read.csv(file = micsname)
micnames <- 1:dim(mics)[1]
mics <- cbind(micnames,mics)
names(mics) <- c("name","x","y")
add <- max(diff(mics$x),diff(mics$y))
trap.xlim <- range(mics$x)+0.5*c(-add,add)
trap.ylim <- range(mics$y)+0.5*c(-add,add)
captures <- read.csv(file = detsname)
traps <- read.traps(data = mics, detector = "proximity")
sstraps <- read.traps(data = mics, detector = "signal")
border <- max(dist(traps))/4
buffer <- border*8 # guess at buffer size - need to experiment with this
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
capt <- make.capthist(captures, traps, fmt = "trapID", noccasions=1)
sscapt <- make.capthist(captures, sstraps, fmt = "trapID", noccasions = 1, cutval = 150)
dists <- distances(traps, mask)

## Setting things up for signal strength analysis.
capt.ss <- capt
for (i in 1:length(captures$ss)){
    capt.ss[captures$ID[i], , captures$trap[i]] = captures$ss[i]
}
