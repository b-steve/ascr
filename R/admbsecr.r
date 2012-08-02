## admbsecr() takes capture history and mask objects from the secr
## package and fits an SECR model using ADMB.
admbsecr <- function(capt, traps, mask, sv = c(2000, 0.9, 10, 5), ssqtoa = NULL,
                     angs = NULL, admbwd = NULL, method = "simple",
                     profpars = NULL){
  require(R2admb)
  require(secr)
  ## Warnings for incurrect input.
  if (length(method) != 1){
    stop("method must be of length 1")
  }
  if (method == "simple" & any(capt != 1 & capt != 0)){
    stop('capt must be binary when using the "simple" method')
  }
  prefix <- paste(method, "secr", sep = "")
  currwd <- getwd()
  ## Moving to ADMB working directory.
  if (!is.null(admbwd)){
    setwd(admbwd)
  }
  ## If NAs are present in capture history object, change to zeros.
  capt[is.na(capt)] <- 0
  ## Extracting no. animals trapped (n) and traps (k) from capture history array.
  ## Only currently works with one capture session.
  n <- dim(capt)[1]
  k <- dim(capt)[3]
  ## Area covered by each mask location.
  A <- attr(mask, "area")
  ## Removing attributes from capt and mask objects as do_admb cannot handle them.
  capt <- matrix(as.vector(capt), nrow = n, ncol = k)
  mask <- as.matrix(mask)
  ## No. of mask locations.
  nm <- nrow(mask)
  ## Distances between traps and mask locations.
  dist <- distances(traps, mask)
  ## Calculating sensible start values.
  if (any(sv == "auto")){
    sv <- autosv(n, nm, A, dist, method)
    print(sv)
  }
  ## Setting up parameters for do_admb.
  if (method == "simple"){
    data <- list(n = n, ntraps = k, nmask = nm, A = A, capt = capt, dist = dist)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3])
    bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000))
  } else if (method == "toa"){
    if (is.null(ssqtoa)){
      ssqtoa <- apply(capt, 1, toa.ssq, dists = dist)
    }
    data <- list(n = n, ntraps = k, nmask = nm, A = A, toacapt = capt,
                 toassq = t(ssqtoa), dist = dist)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], sigmatoa = sv[4])
    bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000), sigmatoa = c(0, 100000))
  } else if (method == "ang"){
    if (is.null(angs)){
      angs <- angles(traps, mask)
    }
    data <- list(n = n, ntraps = k, nmask = nm, A = A, angcapt = capt, ang = angs, dist = dist)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], kappa = sv[4])
    bounds <- list(D = c(0, 100000), g0 = c(0, 1), sigma = c(0, 100000), kappa = c(0, 100000))
  } else {
    stop('method must be either "simple" or "toa"')
  }
  ## Fitting the model.
  if (!is.null(profpars)){
    fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = TRUE,
                   profile = TRUE, profpars = profpars,
                   run.opts = run.control(checkdata = "write", checkparam = "write", clean = TRUE))
  } else {
    fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = TRUE,
                   run.opts = run.control(checkdata = "write", checkparam = "write", clean = TRUE))
  }
  setwd(currwd)
  fit
}

autosv <- function(n, nm, A, dist, method){
  g0 <- 0.9
  sigma <- sqrt(-max(dist)^2/(2*log(1e-100/g0)))
  D <- n/(sum(pdot(mask, traps, 0, list(g0 = g0, sigma = sigma), 1))*A)
  sv <- c(D, g0, sigma)
  if (method == "toa"){
    sigmatoa <- 0.0025
    sv <- c(sv, sigmatoa)
  } else if (method == "ang"){
    kappa <- 10
    sv <- c(sv, kappa)
  }
  sv
}



g0 <- 1
sigma <- 3
t1 <- -1
t2 <- 1

x <- seq(-20, 20, length.out = 500)
d1 <- g0*exp(-(x - t1)^2/(2*sigma^2))
d2 <- g0*exp(-(x - t2)^2/(2*sigma^2))

plot(x, d1, type = "l", col = "blue")
lines(x, d2, col = "red")
lines(x, d1*d2, col = "green")

sh <- numeric(1000)
for (j in 1:1000){
  x <- runif(100000, -20, 20)
  p1 <- g0*exp(-(x - t1)^2/(2*sigma^2))
  p2 <- g0*exp(-(x - t2)^2/(2*sigma^2))
  
  s1 <- numeric(100000)
  s2 <- numeric(100000)
  
  for (i in 1:100000){
    s1[i] <- rbinom(1, 1, p1[i])
    s2[i] <- rbinom(1, 1, p2[i])
  }
  
  n1 <- sum(s1)
  n2 <- sum(s2)
  n3 <- sum(s1*s2)
  
  s[j] <- sigmahat(n1, n2, n3, t)
}

defint <- function(sigma){
  2*sigma*sqrt(pi/2)
}

jointdefint <- function(sigma, t){
  sigma*sqrt(pi)*exp(-t^2/sigma^2)*exp(t^2/sigma^2)^0.75
}

sigmahat <- function(n1, n2, n3, t){
  n <- n1 + n2
  N <- 2*n3
  sqrt(-t^2/(4*log(N*sqrt(2)/n)))
}
