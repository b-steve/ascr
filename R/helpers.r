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
        ID <- ID+1
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

#' Sum of Squares TOA matrix
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
