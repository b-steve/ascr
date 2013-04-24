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
ssmrdsres <- read.table(paste(resfile, "ssmrdsres.txt", sep = ""), header = TRUE)

resmats <- list(hnres, hnfixres, hrres, hrfixres, ssres, logssres, thres, mrdsres, ssmrdsres)

probs <- NULL
for (i in 1:length(resmats)){
  probs <- c(probs, which(resmats[[i]]$maxgrad < -1))
}

hnD <- hnres$D[hnres$maxgrad > -1 & !is.na(hnres$D)]
hnfixD <- hnfixres$D[hnfixres$maxgrad > -1 & !is.na(hnfixres$D)]
hrD <- hrres$D[hrres$maxgrad > -1 & !is.na(hrres$D)]
hrfixD <- hrfixres$D[hrfixres$maxgrad > -1 & !is.na(hrfixres$D)]
ssD <- ssres$D[ssres$maxgrad > -1 & !is.na(ssres$D)]
logssD <- logssres$D[logssres$maxgrad > -1 & !is.na(logssres$D)]
thD <- thres$D[thres$maxgrad > -1 & !is.na(thres$D)]
mrdsD <- mrdsres$D[mrdsres$maxgrad > -1 & !is.na(mrdsres$D)]
ssmrdsD <- ssmrdsres$D[ssmrdsres$maxgrad > -1 & !is.na(ssmrdsres$D)]



dhnD <- density(hnD)
dhnfixD <- density(hnfixD)
dhrD <- density(hrD)
dhrfixD <- density(hrfixD)
dssD <- density(ssD)
dlogssD <- density(logssD)
dthD <- density(thD)
dmrdsD <- density(mrdsD)
dssmrdsD <- density(ssmrdsD)
xs <- c(dhnD$x, dhnfixD$x, dhrD$x, dhrfixD$x, dssD$x,
        dlogssD$x, dthD$x, dmrdsD$x, dssmrdsD$x)
ys <- c(dhnD$y, dhnfixD$y, dhrD$y, dhrfixD$y, dssD$y,
        dlogssD$y, dthD$y, dmrdsD$y, dssmrdsD$y)

##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(dhnD, col = "blue")
#lines(dhnfixD) # These two the same.
lines(dhrD, col = "red")
lines(dhrfixD, col = "purple")
lines(dssD, col = "orange")
lines(dlogssD, col = "darkseagreen")
lines(dthD, col = "brown")
lines(dmrdsD, col = "green")
lines(dssmrdsD, col = "black")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density",
      xlab = expression(hat(D)), ylab = "Density")
legend(x = "topright", legend = c("Half normal", "Hazard rate", "Hazard rate (g0 = 1)",
                         "SS", "SS (log link)", "Threshold", "MRDS with threshold", "MRDS with SS"),
       col = c("blue", "red", "purple", "orange", "darkseagreen", "brown", "green", "black"),
       lty = 1)
##dev.off()

## Calculating bias
hn.bias <- mean(hnD) - D
hnfix.bias <- mean(hnfixD) - D
hr.bias <- mean(hrD) - D
hrfix.bias <- mean(hrfixD) - D
ss.bias <- mean(ssD) - D
logss.bias <- mean(logssD) - D
th.bias <- mean(thD) - D
mrds.bias <- mean(mrdsD) - D
ssmrds.bias <- mean(ssmrdsD) - D

## Calculating variance
hn.var <- var(hnD)
hnfix.var <- var(hnfixD)
hr.var <- var(hrD)
hrfix.var <- var(hrfixD)
ss.var <- var(ssD)
logss.var <- var(logssD)
th.var <- var(thD)
mrds.var <- var(mrdsD)
ssmrds.var <- var(ssmrdsD)

## Calculating mean square error
hn.mse <- hn.bias^2 + hn.var
hnfix.mse <- hnfix.bias^2 + hnfix.var
hr.mse <- hr.bias^2 + hr.var
hrfix.mse <- hrfix.bias^2 + hrfix.var
ss.mse <- ss.bias^2 + ss.var
logss.mse <- logss.bias^2 + logss.var
th.mse <- th.bias^2 + th.var
mrds.mse <- mrds.bias^2 + mrds.var
ssmrds.mse <- ssmrds.bias^2 + ssmrds.var

sort(c(hn = hn.mse, hnfix = hnfix.mse, hr = hr.mse, hrfix = hrfix.mse, ss = ss.mse,
       logss = logss.mse, th = th.mse, mrds = mrds.mse, ssmrds = ssmrds.mse))

