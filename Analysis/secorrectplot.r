## Plots and analysis results from secorrect simulations.

library(devtools)
load_all("~/admbsecr")
library(secr)
library(R2admb)
source("~/admbsecr/Results/secorrect/1/pars.r")

## Determining the fit identifications that were successful.
setwd("~/admbsecr/Results/secorrect/1/fits")
fit.files <- paste("fits", 1:500, "RData", sep = ".")
nfits <- length(fit.files)

Ds <- numeric(500)
Ds.c <- numeric(500)
ses <- numeric(500)
ses.c <- numeric(500)
shs <- numeric(500)
shs.c <- numeric(500)
shs.ses <- numeric(500)
shs.ses.c <- numeric(500)
scs <- numeric(500)
scs.c <- numeric(500)
scs.ses <- numeric(500)
scs.ses.c <- numeric(500)
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
  shs[j] <- out$fit$coefficients[2]
  shs.c[j] <- out$fit.sec$se.correct$coefficients.corrected[2]
  shs.ses[j] <- out$fit$se[2]
  shs.ses.c[j] <- out$fit.sec$se.correct$se.corrected[2]
  scs[j] <- out$fit$coefficients[3]
  scs.c[j] <- out$fit.sec$se.correct$coefficients.corrected[3]
  scs.ses[j] <- out$fit$se[3]
  scs.ses.c[j] <- out$fit.sec$se.correct$se.corrected[3]
  boots.l[[j]] <- boots
  j <- j + 1
}

mean(Ds)
mean(Ds.c)
sd(Ds)
sd(Ds.c)
mean(ses)
mean(ses.c)

mean(shs)
mean(shs.c)
sd(shs)
sd(shs.c)
mean(shs.ses)
mean(shs.ses.c)

mean(scs)
mean(scs.c)
sd(scs)
sd(scs.c)
mean(scs.ses)
mean(scs.ses.c)

cis.wald <- matrix(0, nrow = 500, ncol = 2)
cis.c.norm <- matrix(0, nrow = 500, ncol = 2)
cis.c.perc <- matrix(0, nrow = 500, ncol = 2)
for (i in 1:500){
  cis.wald[i, ] <- Ds[i] + c(-1, 1)*qnorm(0.975)*ses[i]
  cis.c.norm[i, ] <- Ds.c[i] + c(-1, 1)*qnorm(0.975)*ses.c[i]
  cis.c.perc[i, ] <- quantile(boots.l[[i]][, 1], c(0.025, 0.975), na.rm = TRUE)
}
cov.wald <- apply(cis.wald, 1, function(x) x[1] <= 1750 & x[2] >= 1750)
mean(cov.wald)
cov.c.norm <- apply(cis.c.norm, 1, function(x) x[1] <= 1750 & x[2] >= 1750)
mean(cov.c.norm)
cov.c.perc <- apply(cis.c.perc, 1, function(x) x[1] <= 1750 & x[2] >= 1750)
mean(cov.c.perc)
