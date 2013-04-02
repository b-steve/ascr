#' Calculating distances between mask points and traps for SECR models
#'
#' Calculates pairwise distances between traps and mask points.
#'
#' @param traps matrix containing trap coordinates.
#' @param mask matrix containing mask point coordinates.
#' @return A matrix.
#' @export
distances <- function (traps, mask) {
  traps <- as.matrix(traps)
  mask <- as.matrix(mask)
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(mask, 1, d)
  }
  t(apply(traps, 1, onerow))
}

## distcode <- '
##   NumericMatrix TRAPS(traps);
##   NumericMatrix MASK(mask);
##   int K = TRAPS.nrow();
##   int M = MASK.nrow();
##   NumericMatrix DISTANCES(K,M);
##   for (int m=0; m<M; m++) {
##   	for (int k=0; k<K; k++) {
##       DISTANCES(k,m) = pow(pow(TRAPS(k,0)-MASK(m,0),2)+pow(TRAPS(k,1)-MASK(m,1),2),0.5);
##     }
##   }
##   return wrap(DISTANCES);
## '
## distances.cpp <- cxxfunction(signature(traps = "numeric", mask = "numeric"),
##                               body = distcode, plugin = "Rcpp")


#' Calculating angles between mask points and traps for SECR models
#'
#' Calculates angles from each trap to each mask point.
#'
#' @param traps matrix containing trap coordinates.
#' @param mask matrix containing mask coordinates.
#' @return A matrix.
#' @export
angles <- function (traps, mask) {
  traps <- as.matrix(traps)
  mask <- as.matrix(mask)
  onerow <- function (xy) {
    d <- function(xy2) {
      denom=sqrt(sum((xy2-xy)^2))
      if(denom!=0) {
        sintheta=(xy2[1]-xy[1])/denom
        theta=asin(sintheta)
      } else theta=0
      if(xy2[2]<xy[2] & xy2[1]>=xy[1]) theta=theta=pi - theta
      if(xy2[2]<xy[2] & xy2[1]<xy[1]) theta=pi - theta
      if(xy2[2]<xy[2] & xy2[1]==xy[1]) theta=pi
      if(xy2[2]>=xy[2] & xy2[1]<xy[1]) theta=2*pi + theta
      return(theta)
    }
    apply(mask, 1, d)
  }
  t(apply(traps, 1, onerow))
}

## angcode <- '
##   NumericMatrix TRAPS(traps);
## 	NumericMatrix MASK(mask);
## 	int K = TRAPS.nrow();
## 	int M = MASK.nrow();
## 	NumericMatrix ANGLES(K,M);
## 	double X;
## 	double Y;
## 	double pi = 3.14159265359;
## 	for (int k=0; k<K; k++) {
## 		for (int m=0; m<M; m++) {
## 			X = MASK(m,0)-TRAPS(k,0);
## 			Y = MASK(m,1)-TRAPS(k,1);
## 			ANGLES(k,m) = atan(X/Y);
## 			if(Y<0) ANGLES(k,m) += pi;
## 			if(X<0 & Y>=0) ANGLES(k,m) += 2*pi;
## 		}
## 	}
## 	return wrap(ANGLES);
## 	'
## angles.cpp <- cxxfunction(signature(traps = "numeric", mask = "numeric") ,
##                            body = angcode, plugin = "Rcpp")

## Returns capture trap numbers.
trapvec <- function(capthist){
  x <- apply(capthist, 3, function(x) sum(x > 0))
  rep(1:length(x), times = x)
}

## Returns capture animal ID numbers.
animalIDvec <- function(capthist){
  x <- c(apply(capthist, 3, function(x) which(x > 0)), recursive = TRUE)
  names(x) <- NULL
  as.character(x)
}

#' Assigning ID numbers to sounds.
#'
#' Identifies recaptures and assigns ID numbers to sounds recorded for an SECR model.
#'
#' Detected sounds are assumed to come from the same animal if times of arrival at
#' different microphones are closer together than the time it would take for sound to
#' travel between these microphones.
#'
#' @param mics a matrix containing the coordinates of trap locations.
#' @param clicks a data frame containing (at least): (i) \code{$tim$}, the precise
#' time of arrival of the received sound, and (ii) \code{$trap} the trap at which
#' the sound was recorded.
#' @param dt a \code{K} by \code{K} matrix (where \code{K} is the number of traps)
#' containing the time taken for sound to travel between each pair of traps.
#' @return A data frame. Specifically, the \code{clicks} dataframe, now with a new
#' variable, \code{ID}.
#' @author David Borchers, Ben Stevenson
#' @export
make.acoustic.captures <- function(mics, clicks, dt){
  K <- dim(mics)[1]
  captures <- clicks
  ct <- rep(-Inf, K)
  ID <- 1
  ct[clicks$trap[1]] <- clicks$tim[1]
  new <- FALSE
  nclicks <- length(clicks$tim)
  for (i in 2:nclicks){
    if (ct[clicks$trap[i]] > -Inf){
      nd <- length(which(ct > -Inf))
      captures$ID[(i - nd):(i - 1)] <- ID
      ct <- rep(-Inf, K)
      ct[clicks$trap[i]] <- clicks$tim[i]
      ID <- ID + 1
      if(i == nclicks) captures$ID[i] <- ID
    }
    else {
      ct[clicks$trap[i]] <- clicks$tim[i]
      ctset <- which(ct > -Inf)
      dts <- dt[ctset, clicks$trap[i]]
      cts <- -(ct[ctset] - clicks$tim[i])
      if (any((cts - dts) > 0)) new <- TRUE
      if (new) {
        nd <- length(which(ct > -Inf)) - 1
        captures$ID[(i - nd):(i - 1)] <- ID
        ct <- rep(-Inf, K)
        ct[clicks$trap[i]] <- clicks$tim[i]
        ID <- ID + 1
        new <- FALSE
        if (i == nclicks) captures$ID[i] <- ID
      } else if(i == nclicks){
        nd <- length(which(ct > -Inf))
        captures$ID[(i - nd + 1):i] <- ID
      }
    }
  }
  captures
}

#' Sum of squares TOA matrix
#'
#' Calculates ssqtoa matrix for a SECR model with TOA information.
#'
#' @param wit capture history.
#' @param dists distances.
#' @return A matrix.
#' @export
toa.ssq <- function(wit, dists) {
  ssq <- function(x) sum((x - mean(x))^2)
  v <- 330 # speed of sound
  wit.na <- wit
  wit.na[wit == 0] <- NA
  delt <- na.omit(as.vector(wit.na) - dists/v)
  toassq <- apply(delt, 2, ssq)
  toassq
}

#' Simulated signal strength capture history matrix
#'
#' Simulating a signal strength capture history matrix. Signal
#' strength detection function uses a log link function, and
#' thus is different to \code{sim.capthist} where
#' \code{detectfn = 10}.
#'
#' @param traps trap locations.
#' @param popn simulated population.
#' @param detectpars detection function parameters.
#' @param log.link logical, if \code{TRUE} a log link functin is used for
#' expected signal strengths.
#' @export
sim.capthist.ss <- function(traps, popn, detectpars, log.link){
  ssb0 <- detectpars$beta0
  ssb1 <- detectpars$beta1
  sigmass <- detectpars$sdS
  c <- detectpars$cutval
  ntraps <- nrow(traps)
  n <- nrow(popn)
  dists <- distances(as.matrix(popn), as.matrix(traps))
  muss <- ssb0 + ssb1*dists
  if (log.link){
    muss <- exp(muss)
  }
  ss.error <- matrix(rnorm(n*ntraps, 0, sigmass), nrow = n, ncol = ntraps)
  ss <- muss + ss.error
  ss[ss < c] <- 0
  dets <- apply(ss, 1, function(x) !all(x == 0))
  ndet <- sum(dets)
  ss <- ss[dets, ]
  array(ss, dim = c(ndet, 1, ntraps), dimnames = list(which(dets), NULL, NULL))
}

#' Simulated distance error model capture history matrix
#'
#' Simulating a capture history matrix for a distance error model.
#' This function specifically deals with the case where there are
#' two traps in the same location, each of which has different
#' detection function parameters.
#'
#' @param traps trap locations.
#' @param popn simulated population.
#' @param detectpars detection function parameters.
#' @export
sim.capthist.dist <- function(traps, popn, detectpars){
  g01 <- detectpars$g01
  g02 <- detectpars$g02
  sigma1 <- detectpars$sigma1
  sigma2 <- detectpars$sigma2
  alpha <- detectpars$alpha
  dists <- distances(as.matrix(popn), as.matrix(traps))
  N <- nrow(dists)
  ntraps <- ncol(dists)
  probs1 <- g01*exp(-dists[, 1]^2/(2*sigma1^2))
  probs2 <- g02*exp(-dists[, 2]^2/(2*sigma2^2))
  probs <- cbind(probs1, probs2)
  rnos <- matrix(runif(N*ntraps), nrow = N, ncol = ntraps)
  bincapt <- rnos <= probs
  class(bincapt) <- "numeric"
  bincapt <- array(bincapt, dim = c(dim(bincapt)[1], 1, dim(bincapt)[2]))
  dets <- apply(bincapt, 1, function(x) !all(x == 0))
  bincapt <- bincapt[dets, , , drop = FALSE]
  captdists <- dists[dets, ]
  distests <- function(x, alpha){
    betas <- alpha/x
    c(rgamma(1, shape = alpha, rate = betas[1]),
      rgamma(1, shape = alpha, rate = betas[2]))
  }
  estdists <- t(apply(captdists, 1, distests, alpha = alpha))
  estdists <- array(estdists, dim = c(dim(estdists)[1], 1, dim(captdists)[2]))
  capthist.dist <- as.array(estdists*bincapt)
  distfn <- function(x){
    n <- length(x)
    x <- x[x > 0]
    rep(mean(x), n)
  }
  avedists <- t(apply(capthist.dist, 1, distfn))
  capthist.mrds <- array(c(bincapt, avedists),
                         dim = c(dim(capthist.dist), 2))
  list(dist = capthist.dist, mrds = capthist.mrds)
}

#' Fitting MRDS models.
#'
#' Fits a mark-recapture distance sampling (MRDS) model, with different detection
#' function parameters for each trap. This is a temporary function; eventually
#' \code{admbsecr()} will be flexible enough to do this.
#'
#' @param capt capture history array.
#' @param mask mask point locations.
#' @param traps trap locations.
#' @param sv start values.
#' @param admb.dir directory containing mrdstrapcovsecr.tpl.
#' @param clean logical, if \code{TRUE} ADMB files are cleaned after fitting of the model.
#' @param verbose logical, if \code{TRUE} ADMB details, along with error messages, are
#' printed to the R session.
#' @param trace logical, if \code{TRUE} parameter values at each step of the fitting
#' algorithm are printed to the R session.
#' @param detfn the detection function to be used. Either half normal (\code{"hn"}),
#' hazard rate (\code{"hr"}) or threshold (\code{"th"}).
#' @export
mrdstrapcov <- function(capt, mask, traps, sv, admb.dir, clean, verbose, trace, detfn = "hn"){
  setwd(admb.dir)
  n <- dim(capt)[1]
  k <- dim(capt)[3]
  A <- attr(mask, "area")
  nm <- nrow(mask)
  dist <- distances(traps, mask)
  capt <- array(as.vector(capt), dim = c(n, k, dim(capt)[4]))
  data <- list(n = n, ntraps = k, nmask = nm, A = A, capt = capt[, , 1],
               dist = dist, indivdist = capt[, , 2], trace = as.numeric(trace))
  if (detfn == "hn"){
    params <- list(D = sv[1], g01 = sv[2], sigma1 = sv[3],
                   g02 = sv[4], sigma2 = sv[5])
    bounds <- list(D = c(0, 1e8), g01 = c(0, 1), sigma1 = c(0, 1e5),
                   g02 = c(0, 1), sigma2 = c(0, 1e5))

  } else if (detfn == "th"){
    params <- list(D = sv[1], shape1 = sv[2], scale1 = sv[3], shape2 = sv[4],
                   scale2 = sv[5])
    bounds <- list(D = c(0, 1e8), scale1 = c(-10, 0), scale2 = c(-10, 0))
  } else if (detfn == "hr"){
    params <- list(D = sv[1], g01 = sv[2], sigma1 = sv[3], z1 = sv[4],
                   g02 = sv[5], sigma2 = sv[6], z2 = sv[7])
    bounds <- list(D = c(0, 1e8), g01 = c(0, 1), sigma1 = c(0, 1e5),
                   g02 = c(0, 1), sigma2 = c(0, 1e5))
  }
  filename <- paste("mrdstrapcovsecr", detfn, sep = "")
  fit <- do_admb(filename, data = data, params = params, bounds = bounds,
                 verbose = verbose, safe = FALSE,
                 run.opts = run.control(checkdata = "write", checkparam = "write",
                   clean_files = clean))
  fit$data <- data
  fit$traps <- traps
  fit$mask <- as.matrix(mask)
  fit$method <- "mrdstc"
  class(fit) <- c(class(fit), fit$method, "admbsecr")
  fit
}

#' Fitting distance error SECR models.
#'
#' Fits a distance error SECR model, with different detection function parameters for each
#' trap. This is a temporary function; eventually \code{admbsecr()} will be flexible enough
#' to do this.
#'
#' @param capt capture history array.
#' @param mask mask point locations.
#' @param traps trap locations.
#' @param sv start values.
#' @param admb.dir directory containing mrdstrapcovsecr.tpl.
#' @param clean logical, if \code{TRUE} ADMB files are cleaned after fitting of the model.
#' @param verbose logical, if \code{TRUE} ADMB details, along with error messages, are
#' printed to the R session.
#' @param trace logical, if \code{TRUE} parameter values at each step of the fitting
#' algorithm are printed to the R session.
#' @param detfn the detection function to be used. Either half normal (\code{"hn"}),
#' hazard rate (\code{"hr"}) or threshold (\code{"th"}).
#' @export
disttrapcov <- function(capt, mask, traps, sv, admb.dir, clean, verbose, trace, detfn = "hn"){
  setwd(admb.dir)
  n <- dim(capt)[1]
  k <- dim(capt)[3]
  A <- attr(mask, "area")
  nm <- nrow(mask)
  dist <- distances(traps, mask)
  capt <- array(as.vector(capt), dim = c(n, k))
  bincapt <- capt
  bincapt[bincapt != 0] <- 1
  data <- list(n = n, ntraps = k, nmask = nm, A = A, distcapt = capt,
               dist = dist, capt = bincapt, trace = as.numeric(trace))
  if (detfn == "hn"){
    params <- list(D = sv[1], g01 = sv[2], sigma1 = sv[3],
                   g02 = sv[4], sigma2 = sv[5], alpha = sv[6])
    bounds <- list(D = c(0, 1e8), g01 = c(0, 1), sigma1 = c(0, 1e5),
                   g02 = c(0, 1), sigma2 = c(0, 1e5), alpha = c(0, 150))

  } else if (detfn == "th"){
    params <- list(D = sv[1], shape1 = sv[2], scale1 = sv[3], shape2 = sv[4],
                   scale2 = sv[5], alpha = sv[6])
    bounds <- list(D = c(0, 1e8), scale1 = c(-10, 0), scale2 = c(-10, 0),
                   alpha = c(0, 150))
  } else if (detfn == "hr"){
    params <- list(D = sv[1], g01 = sv[2], sigma1 = sv[3], z1 = sv[4],
                   g02 = sv[5], sigma2 = sv[6], z2 = sv[7], alpha = sv[8])
    bounds <- list(D = c(0, 1e8), g01 = c(0, 1), sigma1 = c(0, 1e5),
                   g02 = c(0, 1), sigma2 = c(0, 1e5), alpha = c(0, 150))
  }
  filename <- paste("disttrapcovsecr", detfn, sep = "")
  fit <- do_admb(filename, data = data, params = params, bounds = bounds,
                 verbose = verbose, safe = FALSE,
                 run.opts = run.control(checkdata = "write", checkparam = "write",
                   clean_files = clean))
  fit$data <- data
  fit$traps <- traps
  fit$mask <- as.matrix(mask)
  fit$method <- "disttc"
  class(fit) <- c(class(fit), fit$method, "admbsecr")
  fit
}
