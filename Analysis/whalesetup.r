setwd(dat.dir)
## Get and manipulate data.
fake.traps <- data.frame(name = 1:2, x = c(1e-3, -1e-3), y = c(-1e-3, 1e-3))
real.traps <- data.frame(name = 1:2, x = rep(0, 2), y = rep(0, 2))
whales87.df <- read.csv("ccdat87.csv", header = TRUE)
whales87.df <- whales87.df[whales87.df$sighting != 9999, ]
old.sighting <- whales87.df$sighting
sighting <- numeric(length(old.sighting))
j <- 1
for (i in unique(old.sighting)){
  sighting[old.sighting == i] <- j
  j <- j + 1
}
capthist87.dist <- matrix(0, nrow = max(sighting), ncol = 2)
for (i in unique(sighting)){
  capthist87.dist[i, ] <- whales87.df$r[sighting == i]
}
capthist87.dist[capthist87.dist == 9999] <- 0
capthist87.dist <- capthist87.dist*1000
capthist87.dist <- array(capthist87.dist,
                         dim = c(nrow(capthist87.dist), 1, ncol(capthist87.dist)))
options(warn = -1)
fake.traps <- read.traps(data = fake.traps, detector = "proximity")
real.traps <- read.traps(data = real.traps, detector = "proximity")
## Buffer given by Pr(X <= x) = 0.999 for longest estimated distance.
buffer <- 3000
mask87 <- make.mask(fake.traps, buffer = buffer, type = "trapbuffer")
options(warn = 1)
bincapt87 <- capthist87.dist
bincapt87[bincapt87 > 0] <- 1

whales01.df <- read.csv("whales2001.csv", header = TRUE)
whales01.df <- whales01.df[whales01.df$sighting != -1, ]
old.sighting <- whales01.df$sighting
sighting <- numeric(length(old.sighting))
j <- 1
for (i in unique(old.sighting)){
  sighting[old.sighting == i] <- j
  j <- j + 1
}
capthist01.dist <- matrix(0, nrow = max(sighting), ncol = 2)
for (i in unique(sighting)){
  platord <- whales01.df$plat[sighting == i]
  capthist01.dist[i, ] <- whales01.df$r[sighting == i][platord]
}
capthist01.dist[is.na(capthist01.dist)] <- 0
capthist01.dist <- capthist01.dist*1000
capthist01.dist <- array(capthist01.dist,
                         dim = c(nrow(capthist01.dist), 1, ncol(capthist01.dist)))
options(warn = -1)
## Buffer given by Pr(X <= x) = 0.999 for longest estimated distance.
buffer <- 4500
mask01 <- make.mask(fake.traps, buffer = buffer, type = "trapbuffer")
options(warn = 1)
bincapt01 <- capthist01.dist
bincapt01[bincapt01 > 0] <- 1
distfn <- function(x){
  n <- length(x)
  x <- x[x > 0]
  rep(mean(x), n)
}
dists01 <- t(apply(capthist01.dist, 1, distfn))
capthist.mrds <- array(c(bincapt01, dists01), dim = c(dim(capthist01.dist), 2))
