setwd(dat.dir)
## Get trap locations.

if (datasource == "Res"){
  mics <- read.csv(file = "Site1_estlocs.csv")
  alldat <- read.csv("20120516site1_data_out_of_Access.csv")[501:1000, ]
} else if (datasource == "Original"){
  mics <- read.csv(file = "array1a-top.csv")
  alldat <- read.csv("array1a-top-data.csv")
}

micnames <- 1:dim(mics)[1]
mics <- cbind(micnames, mics)
names(mics) <- c("name", "x", "y")
add <- max(diff(mics$x), diff(mics$y))
trap.xlim <- range(mics$x) + 0.5*c(-add, add)
trap.ylim <- range(mics$y) + 0.5*c(-add, add)

dat <- alldat
## Convert times to milliseconds.
dat$tim <- alldat$startSeconds*1000
## Start at time 0.
dat$tim <- dat$tim - dat$tim[1]
dat$mic <- log2(alldat$channelMap) + 1
dat$ss <- alldat$amplitude
keep <- c("tim","mic","ss")
dat <- dat[, keep]
## Extract time and microphone and arrange in order of increasing time.
ord <- order(dat$tim)
dat <- dat[ord, ]
n <- length(dat$tim)
## Create objects for secr analysis.
## Speed of sound.
sspd <- 330/1000
dloc <- as.matrix(dist(mics))
## Time it takes sound to get between each microphone.
dt <- dloc/sspd
## Make trap object for secr analysis.
traps <- read.traps(data = mics, detector = "proximity")
sstraps <- read.traps(data = mics, detector = "signal")
## Set border for plotting to be 1/4 max distance between traps.
border <- max(dist(traps))/4
if (datasource == "Res"){
  buffer <- border*10
} else if (datasource == "Original"){
  buffer <- border*8
}
## Make capture history object with times of detection
clicks <- data.frame(session = rep(1, n), ID = 1:n, occasion = rep(1, n), trap = dat$mic,
                     tim = dat$tim, ss = dat$ss)
## Reset times to start at 1.
clicks$tim <- clicks$tim - clicks$tim[1] + 1
## Clicks reformatted as input for make.capthist, with duplicate ID rule.
captures <- make.acoustic.captures(traps, clicks, dt)
captures[5:6] <- captures[6:5]
colnames(captures)[5:6] <- colnames(captures)[6:5]
## Make capture history objects.
capt <- make.capthist(captures, traps, fmt = "trapID", noccasions = 1)
sscapt <- make.capthist(captures, sstraps, fmt = "trapID", noccasions = 1, cutval = 150)
captures$ss <- clicks$ss
mask <- make.mask(traps, buffer = buffer, type = "trapbuffer")

## Setting things up for TOA analysis.
capt.toa <- capt
## Put TOA into capture history.
for(i in 1:length(captures$tim)){
  capt.toa[captures$ID[i], , captures$trap[i]] <- captures$tim[i]
}
## Conversion to seconds.
capt.toa <- capt.toa/1000
## Pre-calculate distances from each trap to each grid point in the mask.
dists <- distances(as.matrix(traps), as.matrix(mask))

## Set things up for signal strength analysis.
capt.ss <- capt
for (i in 1:length(captures$ss)){
    capt.ss[captures$ID[i], , captures$trap[i]] = captures$ss[i]
}

## Set things up for a joint signal strength and TOA analysis.
capt.joint <- array(c(capt.ss, capt.toa), dim = c(dim(capt), 2))
