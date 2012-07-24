library("secr")
library("CircStats")
#library("clim.pact")
source("/home/ben/SECR/R/helpers.r")
source("/home/ben/SECR/R/admbsecr.r")
load("/home/ben/SECR/Data/Gibbons/gibbons_data.RData")
gibbons <- read.table(file = "/home/ben/SECR/Data/Gibbons/gibbons.txt", header = TRUE)
npoints <- length(unique(gibbons$point))
ntraps <- 3
traps <- make.grid(nx=ntraps, ny=1, spacing=500, detector="proximity")
rownames(traps) <- c("A","B","C")
buffer <- 5000
spacing <- 50
mask <- make.mask(traps, buffer, spacing, type="trapbuffer") ; head(mask) ; dim(mask)
area.ha <- attr(mask, 'area') * nrow(mask)
mask.dists <- distances(traps, mask)
mask.angs  <- angles(traps, mask)
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

p <- with(start.hn, c(log(D),logit(g0),log(sigma),log(kappa)))
##time1 <- system.time({fit <- nlm(f = secrlikelihood.angs.dk.v1, detectfn = 0, g0.fixed = F,
                                 p = p, capthist = radians, mask = mask, dists = mask.dists,
                                 angs = mask.angs, trace=TRUE)})
##time2 <- system.time({fit2 <- admbsecr(capt = radians, traps = traps, mask = mask,
                                       sv = c(0.1125153, 0.95, 750, 10), angs = mask.angs,
                                       admbwd = "/home/ben/SECR/ADMB", method = "ang")})

