## Functions that automatically generate starting values for parameters.

## Not lifted from the secr package.
## Grabs the average recapture distance, or something.
autosigma2 <- function(capt, traps){
  bincapt <- capt$bincapt
  ave.rc.dist <- function(x){
    trap.ids <- which(x == 1)
    if (length(trap.ids) > 1){
      rc.locs <- traps[trap.ids, ]
      rc.dists <- distances(rc.locs, rc.locs)
      w <- length(rc.dists[rc.dists > 0])
      out <- mean(rc.dists[rc.dists > 0])
    } else {
      out <- NA
      w <- NA
    }
    c(out, w)
  }
  mean.dists <- apply(bincapt, 1, ave.rc.dist)
  mean.dists <- mean.dists[, !is.na(mean.dists[1, ])]
  sum(mean.dists[1, ]*mean.dists[2, ])/sum(mean.dists[2, ])
}

## Write own pdot function.
autoD2 <- function(capt, traps){
  1000
}

autog02 <- function(capt, traps){
  0.95
}

autosigma.toa2 <- function(capt, traps){
  0.0025
}

autokappa2 <- function(capt, traps){
  10
}

autob0.ss2 <- function(capt, traps){
  143
}

autob1.ss2 <- function(capt, traps){
  1
}

autosigma.ss2 <- function(capt, traps){
  5.7
}

autoalpha2 <- function(capt, traps){
  2
}


autoshape2 <- function(capt, traps){
  5
}

autoscale2 <- function(capt, traps){
  10
}
