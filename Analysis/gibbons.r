library("secr")
library("CircStats")
library("inline")
library("Rcpp")

if (.Platform$OS == "unix"){
    source("/home/ben/SECR/R/helpers.r")
    source("/home/ben/SECR/R/admbsecr.r")
    load("/home/ben/SECR/Data/Gibbons/gibbons_data.RData")
    admb.dir <- "/home/ben/SECR/ADMB"
    dat.dir <- "/home/ben/SECR/Data/Gibbons/gibbons.txt"
} else if (.Platform$OS == "windows"){
    source("C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R\\helpers.r")
    source("C:\\Documents and Settings\\Ben\\My Documents\\SECR\\R\\admbsecr.r")
    load("C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Gibbons\\gibbons_data.RData")
    admb.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\ADMB"
    dat.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\SECR\\Data\\Gibbons\\gibbons.txt"
}

gibbons <- read.table(file = dat.dir, header = TRUE)
npoints <- length(unique(gibbons$point))
ntraps <- 3
traps <- make.grid(nx=ntraps, ny=1, spacing=500, detector="proximity")
rownames(traps) <- c("A","B","C")
buffer <- 5000
spacing <- 50
mask <- make.mask(traps, buffer, spacing, type="trapbuffer") ; head(mask) ; dim(mask)
area.ha <- attr(mask, 'area') * nrow(mask)
mask.dists <- distances(traps, mask)
##mask.angs  <- system.time(angles(traps, mask))
mask.angs <- angles.cpp(as.matrix(traps), as.matrix(mask))
ndets <- nrow(gibbons)
ncues <- length(unique(gibbons$group.id))
radians <- array(NA, dim=c(ncues,1,ntraps),
                 dimnames=list(group=unique(gibbons$group.id),
                   occasion=1,post=c("A","B","C")))
for(det in 1:ndets){ # det=1
  cue <- which(dimnames(radians)$group==gibbons$group.id[det])
  trap <- which(dimnames(radians)$post==gibbons$post[det])
  radians[cue,1,trap] <- gibbons$bearing[det]/360*2*pi
}
capthist <- array(1, dim(radians), dimnames(radians))
capthist[which(is.na(radians), arr.ind=T)] <- 0
nrecaps <- sum(rowSums(capthist)>1)
## detection function parameters
start.hn <- list(g0=0.95, sigma=750)
start.hr <- list(g0=0.95, sigma=1000, z=5)
## density
start.hn$D <- ncues/(sum(pdot(mask, traps, 0, start.hn, 1))*attr(mask, 'area'))
start.hr$D <- ncues/(sum(pdot(mask, traps, 1, start.hr, 1))*attr(mask, 'area'))
## angles distribution shape parameter
start.hn$kappa <- start.hr$kappa <- 10 # von Mises
start.hn$rho <- start.hr$rho <- 0.75 # wrapped Cauchy
n <- ncues
S <- 1 
K <- ntraps
A <- attr(mask, "area")
M <- nrow(mask)
capthist <- radians
capthist[capthist > 0] <- 1
capthist[is.na(capthist)] <- 0
hash1 <- which(capthist[,1,]==1, arr.ind=T)-1
hash0 <- which(capthist[,1,]==0, arr.ind=T)-1
##mask.dists <- distances(traps, mask)
mask.dists <- distances.cpp(as.matrix(traps), as.matrix(mask))

p <- with(start.hn, c(log(D),logit(g0),log(sigma),log(kappa)))
##time1 <- system.time({fit <- nlm(f = secrlikelihood.angs.dk.v1, detectfn = 0, g0.fixed = F,
##                                 p = p, capthist = radians, mask = mask, dists = mask.dists,
##                                 angs = mask.angs, trace=TRUE)})

## For same start values:
sv <- c("D" = start.hn$D, "g0" = start.hn$g0, "sigma" = start.hn$sigma,
        "kappa" = start.hn$kappa)
time2 <- system.time({fit2 <- admbsecr(capt = radians, traps = traps, mask = mask,
                                       sv = "auto", angs = mask.angs,
                                       admbwd = admb.dir, method = "ang", autogen = FALSE, trace = TRUE)})

time3 <- system.time({fit3 <- nlm(f = secrlikelihood.cpp, p = p, method = 1, ncues = n,
                                  ntraps = K, npoints = M, radians = radians[, 1, ],
                                  hash1 = hash1, hash0 = hash0, mask_area = A,
                                  mask_dists = mask.dists, mask_angs = mask.angs,
                                 hessian = TRUE)})
