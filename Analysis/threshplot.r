## To read in simulation results from a file.
admbsecr.dir <- "~/admbsecr" ## Point this to admbsecr folder.
resfile <- paste(admbsecr.dir, "Results/thresh/1/", sep = "/")
source(paste(resfile, "pars.r", sep = ""))
hnres <- read.table(paste(resfile, "hnres.txt", sep = ""), header = TRUE)
hnfixres <- read.table(paste(resfile, "hnfixres.txt", sep = ""), header = TRUE)
hrres <- read.table(paste(resfile, "hrres.txt", sep = ""), header = TRUE)
hrfixres <- read.table(paste(resfile, "hrfixres.txt", sep = ""), header = TRUE)
ssres <- read.table(paste(resfile, "ssres.txt", sep = ""), header = TRUE)
logssres <- read.table(paste(resfile, "logssres.txt", sep = ""), header = TRUE)
thres <- read.table(paste(resfile, "thres.txt", sep = ""), header = TRUE)
mrdsres <- read.table(paste(resfile, "mrdsres.txt", sep = ""), header = TRUE)

resmats <- list(hnres, hnfixres, hrres, hrfixres, ssres, logssres, thres, mrdsres)

probs <- NULL
for (i in 1:length(resmats)){
  probs <- c(probs, which(resmats[[i]]$maxgrad < -1))
}

dhnD <- density(hnres$D)
dhnfixD <- density(hnfixres$D)
dhrD <- density(hrres$D[hrres$maxgrad > -1], na.rm = TRUE)
dhrfixD <- density(hrfixres$D[hrfixres$maxgrad > -1], na.rm = TRUE)
dssD <- density(ssres$D[ssres$maxgrad > -1])

plot(dhnD, xlim = c(3000, 8000))
lines(dhnfixD)
lines(dhrD)
lines(dhrfixD)
lines(dssD)
