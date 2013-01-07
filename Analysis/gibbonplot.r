## To read in simulation results from a file.
admbsecr.dir <- "~/admbsecr" ## Point this to admbsecr folder.
resfile <- paste(admbsecr.dir, "Results/gibbons/4/", sep = "/")
source(paste(resfile, "pars.r", sep = ""))
angres <- read.table(paste(resfile, "angres.txt", sep = ""), header = TRUE)
simpleres <- read.table(paste(resfile, "simpleres.txt", sep = ""), header = TRUE)

## Assigning the columns to vectors.
for (i in colnames(simpleres)){
  name <- paste("sim", i, sep = "")
  assign(name, simpleres[, i])
}
for (i in colnames(angres)){
  name <- paste("ang", i, sep = "")
  assign(name, angres[, i])
}

## Two different bandwidth selections.
dsimD <- density(simD)
dangD <- density(angD)
##dsimD <- density(simD, bw = "bcv")
##dangD <- density(angD, bw = "bcv")
xs <- c(dsimD$x, dangD$x)
ys <- c(dsimD$y, dangD$y)

##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(dsimD, col = "blue")
lines(dangD, col = "red")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density", xlab = expression(hat(D)), ylab = "Density")
legend(x = "topright", legend = c("SECR", "SECR + angles"), col = c("blue", "red"), lty = 1)
##dev.off()
