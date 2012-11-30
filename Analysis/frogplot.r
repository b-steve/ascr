## To read in simulation results from a file.
admbsecr.dir <- "~/admbsecr" ## Point this to admbsecr folder.
resfile <- paste(admbsecr.dir, "Results/frogs/1/", sep = "/")
source(paste(resfile, "pars.r", sep = ""))
simpleres <- read.table(paste(resfile, "simpleres.txt", sep = ""), header = TRUE)
toares <- read.table(paste(resfile, "toares.txt", sep = ""), header = TRUE)
ssres <- read.table(paste(resfile, "ssres.txt", sep = ""), header = TRUE)
jointres <- read.table(paste(resfile, "jointres.txt", sep = ""), header = TRUE)

## Assigning the columns to vectors.
for (i in colnames(simpleres)){
  name <- paste("sim", i, sep = "")
  assign(name, simpleres[, i])
}
for (i in colnames(toares)){
  name <- paste("toa", i, sep = "")
  assign(name, toares[, i])
}
for (i in colnames(ssres)){
  name <- paste("ss", i, sep = "")
  assign(name, ssres[, i])
}
for (i in colnames(jointres)){
  name <- paste("joint", i, sep = "")
  assign(name, jointres[, i])
}

## Calculating densities.
dsimD <- density(simD)
dtoaD <- density(toaD)
dssD <- density(ssD)
djointD <- density(jointD)
xs <- c(dsimD$x, dtoaD$x, dssD$x, djointD$x)
ys <- c(dsimD$y, dtoaD$y, dssD$y, djointD$y)



##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(dsimD, col = "blue")
lines(dtoaD, col = "red")
lines(dssD, col = "green")
lines(djointD, col = "purple")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density",
      xlab = expression(hat(D)), ylab = "Density")
legend(x = "topright", legend = c("SECR", "SECR + TOA", "SECR + SS", "SECR + TOA + SS"),
       col = c("blue", "red", "green", "purple"), lty = 1)
##dev.off()
