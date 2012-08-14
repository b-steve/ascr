if (.Platform$OS == "unix"){
    func.dir <- "/home/ben/SECR/R" # change this to wherever secracfuns.r and frogID.r are kept.
    dat.dir <- "/home/ben/SECR/Data/Silvermine" # change this to wherever data are kept.
    admb.dir <- "/home/ben/SECR/ADMB"
} else if (.Platform$OS == "windows") {
    func.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R"
    dat.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Silvermine"
    admb.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\ADMB"
}

library(secr)
setwd(func.dir)
source("admbsecr.r")
source("helpers.r")
source("frogID.r")

setwd(dat.dir)

## Get and manipulate data.
mics <- read.csv(file = "Silvermine.mics1.csv")
micnames <- 1:dim(mics)[1]
mics <- cbind(micnames,mics)
names(mics) <- c("name","x","y")
add <- max(diff(mics$x),diff(mics$y))
trap.xlim <- range(mics$x)+0.5*c(-add,add)
trap.ylim <- range(mics$y)+0.5*c(-add,add)
captures <- read.csv("Silvermine.dets1.csv")
traps <- read.traps(data = mics, detector = "proximity")
border <- max(dist(traps))/4
buffer <- border*8 # guess at buffer size - need to experiment with this
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")
capt <- make.capthist(captures, traps, fmt = "trapID", noccasions=1)
dists <- distances(traps, mask)

## Setting things up for signal strength analysis.
capt.ss <- capt
for (i in 1:length(captures$ss)){
    capt.ss[captures$ID[i], , captures$trap[i]] = captures$ss[i]
}

## Carry out simple SECR analysis
ts1 <- system.time({fit = secr.fit(capt,model=list(D~1, g0~1, sigma~1),
                    mask = mask, verify = FALSE)})
ts2 <- system.time({fit2 = admbsecr(capt, traps = traps, mask, sv = "auto",
                   admbwd = admb.dir, method = "simple")})
