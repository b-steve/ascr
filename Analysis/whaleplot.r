## To read in simulation results from a file.
admbsecr.dir <- "~/admbsecr" ## Point this to admbsecr folder.
resfile <- paste(admbsecr.dir, "Results/whales/2/", sep = "/")
source(paste(resfile, "pars.r", sep = ""))
## Results matrices:
## distres for model with distance error.
## mrdsres for mrds model.
distres <- read.table(paste(resfile, "distres.txt", sep = ""), header = TRUE)
mrdsres <- read.table(paste(resfile, "mrdsres.txt", sep = ""), header = TRUE)

## Assigning the columns to vectors and calculating densities.  i.e.,
## distD and ddistD are values and density for parameter D for the
## distance error model.
for (i in colnames(distres)){
  name <- paste("dist", i, sep = "")
  assign(name, distres[, i])
  dname <- paste("ddist", i, sep = "")
  assign(dname, density(get(name)))
}
for (i in colnames(mrdsres)){
  name <- paste("mrds", i, sep = "")
  assign(name, mrdsres[, i])
  dname <- paste("dmrds", i, sep = "")
  assign(dname, density(get(name)))
}

xs <- c(ddistD$x, dmrdsD$x)
ys <- c(ddistD$y, dmrdsD$y)

## Uncomment pdf() and dev.off() to create pdf of plot.

##pdf(file = paste(resfile, "fig", sep = ""))
plot.new()
plot.window(xlim = range(xs), ylim = c(0, max(ys)))
axis(1)
axis(2, las = 1)
abline(v = D, lty = "dotted")
lines(ddistD, col = "blue")
lines(dmrdsD, col = "red")
abline(h = 0, col = "grey")
box()
title(main = "Simulated sampling distributions of animal density",
      xlab = expression(hat(D)), ylab = "Density")
legend(x = "topright", legend = c("SECR + Distances", "MRDS"),
       col = c("blue", "red"), lty = 1)
##dev.off()
