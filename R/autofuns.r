## Functions that automatically generate starting values for parameters.

## Lifted from the secr package.
autosigma <- function(capthist = NULL, bincapt, traps, mask, sv = NULL, cutoff = NULL,
                      method = NULL){
  obsRPSV <- RPSV.mod(bincapt, traps)
  secr:::naivesigma(obsRPSV, traps, mask, 0, 1)
}

## Lifted from the secr package.
autoD <- function(capthist = NULL, bincapt, traps, mask, sv, cutoff = NULL, method = NULL){
  n <- dim(bincapt)[1]
  A <- attr(mask, "area")
  if (method == "ss" | method == "sstoa"){
    g0 <- 0.95
    sigma <- autosigma(capthist, bincapt, traps, mask, sv, method)
  } else {
    g0 <- sv["g0"]
    sigma <- sv["sigma"]
  }
  n/(sum(pdot(mask, traps, 0,
              list(g0 = g0, sigma = sigma), 1))*A)
}

## Don't think we can estimate this in advance for a single session.
autog0 <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                   cutoff = NULL, method = NULL){
  0.95
}

## Take average std dev across all individuals' TOAs. This will be an
## overestimate as it assumes sound travels instantaneously. Does not
## seem to work well with frog data.
autosigmatoa <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                         cutoff = NULL, method = NULL){
  ## if (method == "sstoa"){
  ##   capthist <- capthist[, , , 2, drop = FALSE]
  ## }
  ## mean(apply(capthist, 1, function(x) sd(x[x != 0])), na.rm = TRUE)
  0.0025
}

## Haven't come up with a good way for this yet.
autokappa <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                      cutoff = NULL, method = NULL){
  10
}

## Assumes signal strength received is equal to signal strength
## produced. Mean of received signal strengths is then converted to a
## mean for a truncated normal.
autossb0 <- function(capthist, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                     cutoff, method){
  if (method == "sstoa"){
    capthist <- capthist[, , , 1, drop = FALSE]
  }
  ss <- capthist[capthist != 0]
  mu <- mean(ss)
  sigma <- sd(ss)
  alpha <- (cutoff - mu)/sigma
  lambda <- dnorm(alpha)/(1 - pnorm(alpha))
  mu + sigma*lambda
}

## Assumes signal strength received is almost equal to signal strength
## produced. Given a small negative value to get away from the
## paramter bound.
autossb1 <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                     cutoff = NULL, method = NULL){
  -0.1
}

## Assumes signal strength received is equal to signal strength
## produced. Std dev of received signal strengths is then converted to
## a std dev for a truncated normal.
autosigmass <- function(capthist, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                        cutoff, method = NULL){
  if (method == "sstoa"){
    capthist <- capthist[, , , 1, drop = FALSE]
  }
  ss <- capthist[capthist != 0]
  mu <- mean(ss)
  sigma <- sd(ss)
  alpha <- (cutoff - mu)/sigma
  lambda <- dnorm(alpha)/(1 - pnorm(alpha))
  delta <- lambda*(lambda - alpha)
  sqrt(sigma^2*(1 - delta))
}

autoalpha <- function(capthist, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                       cutoff = NULL, method = NULL){
  2
}

## Following are helper functions for naive sigma estimation.
RPSV.mod <- function(capthist, traps){
  w <- split(trapvec(capthist), animalIDvec(capthist))
  temp <- lapply(w, RPSVx, traps = traps)
  temp <- matrix(unlist(temp), nrow = 3)
  temp <- apply(temp, 1, sum, na.rm = TRUE)
  temp <- sqrt((temp[2] + temp[3])/(temp[1] - 1))
  attr(temp, "names") <- NULL
  temp
}

RPSVx <- function(cx, traps) {
  cx <- abs(cx)
  x <- traps$x[cx]
  y <- traps$y[cx]
  n <- length(x)
  c(n = n - 1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) -
    (sum(y))^2/n)
}
RPSVxy <- function(xy) {
  x <- xy$x
  y <- xy$y
  n <- length(x)
  c(n = n - 1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) -
    (sum(y))^2/n)
}
