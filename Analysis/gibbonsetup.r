if (.Platform$OS == "unix"){
    func.dir <- "/home/ben/admbsecr/R" # change this to wherever R functions are kept.
    admb.dir <- "/home/ben/admbsecr/ADMB"
    dat.dir <- "/home/ben/admbsecr/Data/Gibbons/gibbons.txt"
} else if (.Platform$OS == "windows"){
    func.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\admbsecr\\R"
    admb.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\admbsecr\\ADMB"
    dat.dir <- "C:\\Documents and Settings\\Ben\\My Documents\\admbsecr\\Data\\Gibbons\\gibbons.txt"
}

library("secr")
library("CircStats")
library("inline")
library("Rcpp")

## Set working directory to that with the functions.
setwd(func.dir)
## Get SECR functions.
source("admbsecr.r")
source("autofuns.r")
source("helpers.r")
source("lhoodfuns.r")
source("tplmake.r")

gibbons <- read.table(file = dat.dir, header = TRUE)
npoints <- length(unique(gibbons$point))
ntraps <- 3
traps <- make.grid(nx = ntraps, ny = 1, spacing = 500, detector = "proximity")
rownames(traps) <- c("A", "B", "C")
buffer <- 5000
spacing <- 50
mask <- make.mask(traps, buffer, spacing, type = "trapbuffer")
mask.dists <- distances.cpp(as.matrix(traps), as.matrix(mask))
mask.angs <- angles.cpp(as.matrix(traps), as.matrix(mask))
ndets <- nrow(gibbons)
ncues <- length(unique(gibbons$group.id))
radians <- array(NA, dim = c(ncues, 1, ntraps),
                 dimnames = list(group = unique(gibbons$group.id),
                   occasion = 1, post = c("A", "B", "C")))
for(det in 1:ndets){
  cue <- which(dimnames(radians)$group == gibbons$group.id[det])
  trap <- which(dimnames(radians)$post == gibbons$post[det])
  radians[cue, 1, trap] <- gibbons$bearing[det]/360*2*pi
}
n <- ncues
S <- 1 
K <- ntraps
A <- attr(mask, "area")
M <- nrow(mask)
capthist <- radians
capthist[capthist > 0] <- 1
capthist[is.na(capthist)] <- 0
hash1 <- which(capthist[, 1, ] == 1, arr.ind = TRUE) - 1
hash0 <- which(capthist[, 1, ] == 0, arr.ind = TRUE) - 1


