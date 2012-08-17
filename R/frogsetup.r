if (.Platform$OS == "unix"){
    func.dir <- "/home/ben/SECR/R" # change this to wherever secracfuns.r and frogID.r are kept.
    dat.dir <- "/home/ben/SECR/Data/Frogs" # change this to wherever data are kept.
    admb.dir <- "/home/ben/SECR/ADMB"
} else if (.Platform$OS == "windows") {
    func.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R"
    dat.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Frogs"
    admb.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\ADMB"
}

# get library
library(secr)
library(Rcpp)
library(inline)

# get frog-specific SECR functions
setwd(func.dir) # set working directory to that with the functions
source("admbsecr.r")
source("helpers.r")
source("frogID.r")

# Get and manipulate data
# -----------------------
setwd(dat.dir) # set working directory to that with the data
# get trap locations:
mics=read.csv(file="array1a-top.csv")  # output from calclocs(): must be in increasing order of microphone number
micnames=1:dim(mics)[1]
mics=cbind(micnames,mics)
names(mics)=c("name","x","y")
add=max(diff(mics$x),diff(mics$y))
trap.xlim=range(mics$x)+0.5*c(-add,add)
trap.ylim=range(mics$y)+0.5*c(-add,add)
alldat=read.csv("array1a-top-data.csv") # get detection data (from PAMguard)
dat=alldat
#dat$tim=dat$UTCMilliseconds
dat$tim=alldat$startSeconds*1000 # times in milliseconds
dat$tim=dat$tim-dat$tim[1] # start at time zero
dat$mic=log2(alldat$channelMap)+1
dat$ss=alldat$amplitude
keep=c("tim","mic","ss")
dat=dat[,keep]
# extract time and microphone and arrange in order of increasing time:
ord=order(dat$tim)
dat=dat[ord,]
n=length(dat$tim)
# Create objects for secr analysis (initially without time-of-arrival (TOA or TDOA D being "difference") data)
sspd=330/1000 # speed of sound in metres per millisecond
dloc=as.matrix(dist(mics))
dt=dloc/sspd # time it takes sound to get between each microphone (in milliseconds)
# make trap object for secr analysis:
traps=read.traps(data=mics,detector="proximity")
border=max(dist(traps))/4 # set border for plotting to be 1/4 max dist between traps
# make capture history object with times of detection
clicks=data.frame(session=rep(1,n),ID=1:n,occasion=rep(1,n),trap=dat$mic,tim=dat$tim,ss=dat$ss)
#clicks=clicks[-1,] # remove 1st record which is zero, as next record is >11,000
clicks$tim=clicks$tim-clicks$tim[1]+1 # reset times so start at 1
captures=make.frog.captures(traps,clicks,dt) # clicks reformatted as input for make.capthis, with dup ID rule
capt=make.capthist(captures,traps,fmt="trapID",noccasions=1) # make capture history object
captures$ss=clicks$ss
buffer=border*8 # guess at buffer size - need to experiment with this
mask=make.mask(traps,buffer=buffer,type="trapbuffer")

## Setting things up for TOA analysis.
capt.toa=capt # make copy of capture object to which times will be added
for(i in 1:length(captures$tim)) capt.toa[captures$ID[i],,captures$trap[i]]=captures$tim[i] # put TOA into capture history
capt.toa=capt.toa/1000 # convert from milliseconds to seconds
dists=distances(traps,mask) # pre-calculate distances from each trap to each grid point in mask

## Setting things up for signal strength analysis.
capt.ss <- capt
for (i in 1:length(captures$ss)){
    capt.ss[captures$ID[i], , captures$trap[i]] = captures$ss[i]
}


## Carry out simple SECR analysis
##ts1 <- system.time({fit = secr.fit(capt,model=list(D~1, g0~1, sigma~1),
##                    mask = mask, verify = FALSE)})
##ts2 <- system.time({fit2 = admbsecr(capt, traps = traps, mask, sv = "auto",
##                      admbwd = NULL, method = "simple",
##                      memory = 1500000, autogen = TRUE)})

## Carry out TOA analysis
##ssqtoa <- apply(capt.toa,1,toa.ssq,dists=dists) # create ssq matrix in advance
##sigma.toa <- 0.0025 # starting value for toa measurement error std. err.
##start.beta <- c(fit$fit$estimate,log(sigma.toa)) # add sigma.toa to parameters to be estimated
##ttoa1 <- system.time({toafit = nlm(f = secrlikelihood.toa1, p = start.beta, capthist=capt.toa,
##                      mask = mask, dists = dists, ssqtoa = ssqtoa, trace = TRUE)})
##start.beta2 <- c(coef(fit2), sigma.toa)
##ttoa2 <- system.time({toafit2 = admbsecr(capt = capt.toa, traps = traps, mask = mask,
##                      sv = "auto", ssqtoa = ssqtoa, admbwd = admb.dir, method = "toa")})

## Carry out signal strength analysis
##tss <- system.time({ssfit = admbsecr(capt = capt.ss, traps = traps, mask = mask, sv = "auto",
##                    admbwd = admb.dir, method = "ss")})

