## Letting R know where everything is.
admbsecr.dir <- "~/admbsecr" # Point this to the admbsecr file.
if (.Platform$OS == "unix"){
  sep <- "/"
} else if (.Platform$OS == "windows") {
  sep <- "\\"
}
admb.dir <- paste(admbsecr.dir, "ADMB", sep = sep)
work.dir <- paste(admbsecr.dir, "Analysis", sep = sep)
func.dir <- paste(admbsecr.dir, "R", sep = sep)
dat.dir <- paste(admbsecr.dir, "Data", sep = sep)

## Running setup code.
setwd(work.dir)
source("gibbonsetup.r")

## Carrying out angle analysis with nlm() (likelihood function in C++).
p <- c(log(0.1125153), logit(0.95), log(750), log(10))
anglefit1 <- nlm(f = secrlikelihood.cpp, p = p, method = 1, ncues = n,
                 ntraps = K, npoints = M, radians = radians[, 1, ],
                 hash1 = hash1, hash0 = hash0, mask_area = A,
                 mask_dists = mask.dists, mask_angs = mask.angs,
                 hessian = TRUE)

## Carrying out angle analysis with admbsecr(). Start values same as
## above provided, need to change to sv = sv.
sv <- c("D" = exp(p[1]), "g0" = invlogit(p[2]), "sigma" = exp(p[3]), "kappa" = exp(p[4]))
anglefit2 <- admbsecr(capt = radians, traps = traps, mask = mask,
                      sv = "auto", admbwd = admb.dir, method = "ang")
