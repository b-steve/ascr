pref <- c("Silvermine")
no <- 1

detsname <- paste(pref, paste("dets", no, sep = ""), "csv", sep = ".")
micsname <- paste(pref, paste("mics", no, sep = ""), "csv", sep = ".")

if (.Platform$OS == "unix"){
    func.dir <- "/home/ben/SECR/R" # change this to wherever secracfuns.r and frogID.r are kept.
    dat.dir <- "/home/ben/SECR/Data/Acoustic" # change this to wherever data are kept.
    admb.dir <- "/home/ben/SECR/ADMB"
} else if (.Platform$OS == "windows") {
    func.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R"
    dat.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Acoustic"
    admb.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\ADMB"
}

library(secr)
library(Rcpp)
library(inline)
setwd(func.dir)
source("admbsecr.r")
source("helpers.r")
source("frogID.r")


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

## Carry out simple SECR analysis
##ts1 <- system.time({fit = secr.fit(capt,model=list(D~1, g0~1, sigma~1),
##                    mask = mask, verify = FALSE)})
##ts2 <- system.time({fit2 = admbsecr(capt, traps = traps, mask = mask, sv = "auto",
##                   admbwd = admb.dir, method = "simple", autogen = FALSE, trace = TRUE)})

## Fitting signal strength analysis with nlm().
##startval <- c(log(13000), 170, 0, log(8))
##ssfit1 <- nlm(f = secrlikelihood.ss, p = startval, capthist = capt.ss, mask = mask,
##              dists = dists, cutoff = 150, trace = TRUE)
##ests <- c(exp(ssfit2$estimate[1]), ssfit2$estimate[2:3], exp(ssfit2$estimate[4]))
## Fitting signal strength analysis with secr.fit().
##tss1 <- system.time({ssfit = secr.fit(sscapt,model=list(D~1, g0~1, sigma~1),
##                       detectfn = 10, mask = mask, verify = FALSE, steptol = 1e-4)})
## Fitting signal strength analysis with admbsecr().
##tss2 <- system.time({ssfit2 <- admbsecr(capt.ss, traps = traps, mask, sv = "auto",
##                                        cutoff = 150, admbwd = admb.dir,
##                                        method = "ss", clean = TRUE, trace = TRUE, 
##                                        autogen = FALSE)})
