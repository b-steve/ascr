setwd(dat.dir)
## Get and manipulate data.
fake.traps <- data.frame(name = 1:2, x = c(1e-3, -1e-3), y = c(-1e-3, 1e-3))
real.traps <- data.frame(name = 1:2, x = rep(0, 2), y = rep(0, 2))
whales.df <- read.csv("ccdat87.csv", header = TRUE)
whales.df <- whales.df[whales.df$sighting != 9999, ]
old.sighting <- whales.df$sighting
sighting <- numeric(length(old.sighting))
j <- 1
for (i in unique(old.sighting)){
  sighting[old.sighting == i] <- j
  j <- j + 1
}
capthist.dist <- matrix(0, nrow = max(sighting), ncol = 2)
for (i in unique(sighting)){
  capthist.dist[i, ] <- whales.df$r[sighting == i]
}
capthist.dist[capthist.dist == 9999] <- 0
capthist.dist <- capthist.dist*1000
capthist.dist <- array(capthist.dist, dim = c(nrow(capthist.dist), 1, ncol(capthist.dist)))
options(warn = -1)
fake.traps <- read.traps(data = fake.traps, detector = "proximity")
real.traps <- read.traps(data = real.traps, detector = "proximity")
## Buffer given by Pr(X <= x) = 0.999 for longest estimated distance.
buffer <- 3000
mask <- make.mask(fake.traps, buffer = buffer, type = "trapbuffer")
options(warn = 1)
bincapt <- capthist.dist
bincapt[bincapt > 0] <- 1

