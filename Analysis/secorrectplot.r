## Plots and analysis results from secorrect simulations.

library(devtools)
load_all("~/admbsecr")
library(secr)
library(R2admb)
source("~/admbsecr/Results/secorrect/1/pars.r")

## Determining the fit identifications that were successful.
setwd("~/admbsecr/Results/secorrect/1/fits")
fit.files <- list.files()

nfits <- length(fit.files)

Ds <- numeric(500)
Ds.c <- numeric(500)
ses <- numeric(500)
ses.c <- numeric(500)
boots.l <- list()

j <- 1
for (i in fit.files){
  load(i)
  boots <- out$fit.sec$se.correct$boots
  coefs <- out$fit$coefficients
  bias <- apply(boots, 2, mean, na.rm = TRUE) - coefs
  out$fit.sec$se.correct$coefficients.corrected <- coefs - bias
  out$fit.sec$se.correct$se.corrected <- apply(boots, 2, sd, na.rm = TRUE)
  Ds[j] <- out$fit$coefficients[1]
  Ds.c[j] <- out$fit.sec$se.correct$coefficients.corrected[1]
  ses[j] <- out$fit$se[1]
  ses.c[j] <- out$fit.sec$se.correct$se.corrected[1]
  boots.l[[j]] <- boots
  j <- j + 1
}

mean(Ds)
mean(Ds.c)
mean(ses, na.rm = TRUE)
mean(ses.c)

cis.norm <- matrix(0, nrow = 500, ncol = 2)
cis.perc <- matrix(0, nrow = 500, ncol = 2)
for (i in 1:500){
  cis.norm[i, ] <- Ds.c[i] + c(-1, 1)*qnorm(0.975)*ses.c[i]
  cis.perc[i, ] <- quantile(boots.l[[i]][, 1], c(0.025, 0.975), na.rm = TRUE)
}
cov.norm <- apply(cis.norm, 1, function(x) x[1] <= 1750 & x[2] >= 1750)
mean(cov.norm)
cov.perc <- apply(cis.perc, 1, function(x) x[1] <= 1750 & x[2] >= 1750)
mean(cov.perc)
