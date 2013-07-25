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
Ds <- NULL
Ds.c <- NULL
ses <- NULL
ses.c <- NULL
work.ids <- NULL
j <- 1
for (i in fit.files){
  load(i)
  if (class(out$fit)[1] != "try-error"){
    work.ids <- c(work.ids, i)
    boots <- out$fit.sec$se.correct$boots
    coefs <- out$fit$coefficients
    bias <- apply(boots, 2, mean, na.rm = TRUE) - coefs
    out$fit.sec$se.correct$coefficients.corrected <- coefs - bias
    out$fit.sec$se.correct$se.corrected <- apply(boots, 2, sd, na.rm = TRUE)
    Ds <- c(Ds, out$fit$coefficients[1])
    Ds.c <- c(Ds.c, out$fit.sec$se.correct$coefficients.corrected[1])
    ses <- c(ses, out$fit$se[1])
    ses.c <- c(ses.c, out$fit.sec$se.correct$se.corrected[1])
    j <- j + 1
  } else {
    file.remove(i)
  }
}

mean(Ds)
mean(Ds.c)
mean(ses, na.rm = TRUE)
mean(ses.c)

Ds.c <- Ds.c[!is.na(Ds.c)]

boots <- out$fit.sec$se.correct$boots
coefs <- out$fit$coefficients
bias <- apply(boots, 2, mean, na.rm = TRUE) - coefs
out.fitsec$secorrect$coefficients.corrected <- coefs - bias
out$fitsec$secorrect$se.corrected <- apply(boots, 2, sd, na.rm = TRUE)
