distances <- function (X, Y) {
  ## X and Y are 2-column matrices of coordinates
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}

angles <- function (X, Y) {
#-------------------------------------------------------------------------------
# X and Y are 2-column matrices of coordinates
# Returns angles (0,360] from points in X to points in Y
#         in matrix of dimension nrow(X) x nrow(Y)
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
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}


secrlikelihood <- function (beta, capthist, mask, dist = NULL, trace = FALSE) {
  ## Compute negative log likelihood for halfnormal proximity model
  ## Murray Efford 2011-02-27

  ## Inputs
  ##     beta  -  parameter vector on link scale (D, g0, sigma)
  ##     capthist - capthist object (n x S x K binary array)
  ##     mask - mask object (M x 2 matrix of x-y coords, with attribute 'area')
  ##     dist - K x M matrix of distances between each detector and each mask point
  ##     trace - logical TRUE for output of each evaluation

  ## where K = number of detectors, M = number of mask points

  if (!all(capthist %in% c(0,1)))
    stop ('secrlikelihood requires binary data')

  ## 'real' parameter values
  D <- exp(beta[1])
  ## if D is modelled, expand here to a vector of length
  ## equal to number of rows in mask
  g0 <- invlogit(beta[2])
  sigma <- exp(beta[3])

  n <- nrow(capthist)        ## number observed
  S <- ncol(capthist)        ## number of occasions
  A <- attr( mask, 'area')   ## area of one cell

  ## distances if not passed as argument
  if (is.null(dist)) {
    traps <- traps(capthist)
    dist <- distances(traps, mask)
  }

  ## precompute probabilities for each detector and each mask point
  ## gk is K x M matrix
  gk <- g0 * exp(-dist^2 / 2 / sigma^2)
  log.gk <- log(gk)
  log.gk1 <- log(1-gk)

  ## probability of being caught at least once if at mask site m
  p.m <- 1 - apply(1-gk, 2, prod) ^ S  ## vector of length M
  sumDp <- sum(p.m * D)                ## scalar

  L1 <- sum ( apply(capthist, 1, Dprwi) )


  L2 <- - n * log(sumDp)
  L3 <- dpois (n, sumDp * A, log = TRUE)
  LL <- L1+L2+L3
  if (trace) {
    cat (D, g0, sigma, LL, '\n')
    flush.console()
  }
  if (!is.finite(LL))
    1e10
  else
    -LL   ## return negative log likelihood
}

Dprwi <- function (wi) {
  ## wi is S x K, log.gk is K x M, prwi.s is S x M
  prwi.s <- wi %*% log.gk + (1-wi) %*% log.gk1
  ## sum log(p) over occasions, result is M-vector
  prwi.s <- apply(prwi.s,2,sum)
  log (sum(D * exp(prwi.s)))
}

toa.ssq=function(wit,dists) {
  #-------------------------------------------------------------------------------
  # returns sum over detectors on which animal was detected, of squared
  # (t0hat-mean(t0hat)), where t0hat_m=toa_m-dists_m/v, # toa_m is time of arrival
  # at detector m, dists_m is distance from detector m to each mask point.
  #
  # Returns vecor of length=number of mask points
  #-------------------------------------------------------------------------------
  ssq=function(x) sum((x-mean(x))^2)
  v=330 # speed of sound
  wit.na=wit; wit.na[wit==0]=NA # mark those with no capture
  delt=na.omit(as.vector(wit.na)-dists/v) # need vector else get non-conformable arrays
  toassq=apply(delt,2,ssq)
  return(toassq)
}

secrlikelihood.toa1 <- function (beta, capthist, mask, dists=NULL, ssqtoa=NULL, trace=FALSE) {
  ## Compute negative log likelihood for halfnormal proximity model with TOA data
  ## Murray Efford 2011-02-27
  ## DLB added TOA 2011-10-15;
  ##     update for efficiency (added ssqtoa): 07/11/11

  ## Limitations -
  ##     halfnormal detection function
  ##     'proximity' detector type
  ##     no deaths
  ##     full likelihood
  ##     no competing risk model
  ##     one session
  ##     no groups, covariates, time variation or trap response
  ##     all detectors used
  ##     link functions D = log, g0 = logit, sigma = log

  ## Inputs
  ##     beta  -  parameter vector on link scale (D, g0, sigma, sigma.toa)
  ##     capthist - capthist object with TOAs insteand of 1s (n x S x K binary array)
  ##     mask - mask object (M x 2 matrix of x-y coords, with attribute 'area')
  ##     dists - K x M matrix of distances between each detector and each mask point
  ##     ssqtoa - M x n matrix: sum over detectors for each detection, of (toa-transmission time)^2
  ##     trace - logical TRUE for output of each evaluation

  ## where K = number of detectors, M = number of mask points

  # DLB removed this check when put TOA into capture histories
  #    if (!all(capthist %in% c(0,1)))
  #        stop ('secrlikelihood requires binary data')

  ## 'real' parameter values
  D <- exp(beta[1])
  ## if D is modelled, expand here to a vector of length
  ## equal to number of rows in mask
  g0 <- invlogit(beta[2])
  sigma <- exp(beta[3])
  sigma.toa=exp(beta[4]) # std. dev. of arrival time measurement error

  n <- nrow(capthist)        ## number observed
  S <- ncol(capthist)        ## number of occasions
  A <- attr( mask, 'area')   ## area of one cell

  ## distances if not passed as argument
  if (is.null(dists)) {
    traps <- traps(capthist)
    dists <- distances(traps, mask)
  }

  ## toassq if not passed as argument
  if (is.null(ssqtoa)) {
    ssqtoa <- apply(capthist,1,toa.ssq,dists=dists)
  }

  ## precompute probabilities for each detector and each mask point
  ## gk is K x M matrix
  gk <- g0 * exp(-dists^2 / 2 / sigma^2)
  log.gk <- log(gk)
  log.gk1 <- log(1-gk)

  ## probability of being caught at least once if at mask site m
  p.m <- 1 - apply(1-gk, 2, prod) ^ S  ## vector of length M
  sumDp <- sum(p.m * D)                ## scalar

  log.Dprwi <- function (wit) {
    wi=(wit>0)*1 # wit has detection times and zeros for non-detection
    ## wi is S x K, log.gk is K x M, prwi.s is S x M
    prwi.s <- wi %*% log.gk + (1-wi) %*% log.gk1
    ## sum log(p) over occasions, result is M-vector
    prwi.s <- apply(prwi.s,2,sum) # log(Pr(wis|X)) = log(product over occasions of Pr(wis|X)): log(Eq (4) of Efford, Borchers, Byrom)
    #        if(sum(wi)>1) log.ft=logft(wit,dists,sigma.toa) # log-likelihood for this animal's TOA at each X
    #        else log.ft=0
    #        nc=sum(wi) # number captures
    #        if(nc>1) log.ft=logft(wit,dists,sigma.toa) # log-likelihood for this animal's TOA at each X
    #        else log.ft=0 # so don't add anything for TOA likelihood component for this animal
    #        log (sum(D * exp(prwi.s+log.ft))) # log("integral" of D*Pr(wi|X)*f(toa) over X))
    #        log (D * exp(prwi.s)) # log(D*Pr(wi|X)) at each X)
    log(D) + prwi.s # log(D*Pr(wi|X)) at each X): S x M matrix
  }

  L1.X=apply(capthist,1,log.Dprwi)
  L1.X.toa=log.ftoa(capthist,ssqtoa,sigma.toa)
  #    L1 <- sum (exp(L1.X+L1.X.toa))
  L1 <- sum(log(apply(exp(L1.X+L1.X.toa),2,sum)))
  L2 <- - n * log(sumDp)
  L3 <- dpois (n, sumDp * A, log = TRUE)
  LL <- L1+L2+L3
  if (trace) {
    cat (D, g0, sigma, sigma.toa, LL, '\n')
    flush.console()
  }
  if (!is.finite(LL))
    1e10
  else
    -LL   ## return negative log likelihood
}

log.ftoa=function(capthist,ssqtoa,sigma.toa) {
  num.nonzeros=function(x) sum(x>0)
  logftoa=ssqtoa*0 # initialize
  M=apply(capthist,1,num.nonzeros)
  M2plus=which(M>1)
  madd=matrix(((1-M[M2plus])*log(sigma.toa)),nrow=dim(logftoa)[1],ncol=length(M2plus),byrow=TRUE)
  logftoa[,M2plus]=ssqtoa[,M2plus]/(-2*sigma.toa^2) + madd
  #  logftoa[,M2plus]=t(ssqtoa[,M2plus]/(-2*sigma.toa^2)) + t((1-M[M2plus])*log(sigma.toa)) # omitting terms without parameters
  return(logftoa)
}

secrlikelihood.angs <- function (beta, capthist, mask, dists=NULL, angs=NULL, trace=FALSE) {
## Compute negative log likelihood for halfnormal proximity model with TOA data
## Murray Efford 2011-02-27
## DLB updated secrlikelihood.toa1 to deal with angles instead of TOA 09/11/11

## Limitations -
##     halfnormal detection function
##     'proximity' detector type
##     no deaths
##     full likelihood
##     no competing risk model
##     one session
##     no groups, covariates, time variation or trap response
##     all detectors used
##     link functions D = log, g0 = logit, sigma = log

## Inputs
##     beta  -  parameter vector on link scale (D, g0, sigma, kappa)
##     capthist - capthist object with angeles  (radians) instead of 1s (n x S x K array)
##     mask - mask object (M x 2 matrix of x-y coords, with attribute 'area')
##     dists - K x M matrix of distances between each detector and each mask point
##     angles - K x M matrix: angles from detectors to each X in mask
##     trace - logical TRUE for output of each evaluation

## where K = number of detectors, M = number of mask points

    ## 'real' parameter values
    D <- exp(beta[1])
    ## if D is modelled, expand here to a vector of length
    ## equal to number of rows in mask
    g0 <- invlogit(beta[2])
    sigma <- exp(beta[3])
    kappa=exp(beta[4]) # std. dev. of arrival time measurement error

    n <- nrow(capthist)        ## number observed
    S <- ncol(capthist)        ## number of occasions
    A <- attr( mask, 'area')   ## area of one cell

    ## distances if not passed as argument
    if (is.null(dists)) {
        traps <- traps(capthist)
        dists <- distances(traps, mask)
    }

    ## angs if not passed as argument
    if (is.null(angs)) {
        angs <- angles(traps(capthist),mask)
    }

    ## precompute probabilities for each detector and each mask point
    ## gk is K x M matrix
    gk <- g0 * exp(-dists^2 / 2 / sigma^2)
    log.gk <- log(gk)
    log.gk1 <- log(1-gk)

    ## probability of being caught at least once if at mask site m
    p.m <- 1 - apply(1-gk, 2, prod) ^ S  ## vector of length M
    sumDp <- sum(p.m * D)                ## scalar

    L1.X=apply(capthist,1,log.Dprwi.ang)
    L1.X.a=apply(capthist,1,log.vmCH,mu=angs,kappa=kappa)
    L1 <- sum(log(apply(exp(L1.X+L1.X.a),2,sum)))
    L2 <- - n * log(sumDp)
    L3 <- dpois (n, sumDp * A, log = TRUE)
    LL <- L1+L2+L3
    if (trace) {
        cat (D, g0, sigma, kappa, LL, '\n')
        flush.console()
    }
    if (!is.finite(LL))
        1e10
    else
        -LL   ## return negative log likelihood
}

## Function to calculate log[D*Bern(wi)])=log[D]+log[Bern(wi)] for capture i for each maskpoint
## Works with angle data as missing values are coded with NA in this case.
log.Dprwi.ang <- function (wit) {
  wi <- (!is.na(wit))*1
  prwi.s <- wi %*% log.gk + (1-wi) %*% log.gk1
  prwi.s <-apply(prwi.s,2,sum)
  log(D) + prwi.s
}

log.vmCH <- function(x,mu,kappa) {
  ## returns log-likelihood for one capture history over all cols in mu - i.e. it sums over the detectors (all k)
  ## uses library(CircStats) for von Mises density
  dvm1=function(mu,x,kappa) {dvm(x,mu,k=kappa)}
  dvm.bycol=function(mu,x,kappa) apply(as.matrix(mu,nrow=1),2,dvm1,x=x,kappa=kappa)
  keep=which(!is.na(x)) ; keep
  ##  keep=(x>0)
  logvmCH=log(dvm.bycol(mu[keep,],x[keep],kappa))
  if(dim(logvmCH)[1]==dim(mu)[2] & dim(logvmCH)[2]==1) logvmCH=t(logvmCH) # undo annoying feature of R: turns 1 x K matrices into K x 1 matrices!
  return(apply(logvmCH,2,sum))
}

## Automatically generates starting value for sigma.



RPSV.mod <- function(capthist, traps){
    w <- split(trapvec(capthist), animalIDvec(capthist))
    temp <- lapply(w, RPSVx)
    temp <- matrix(unlist(temp), nrow = 3)
    temp <- apply(temp, 1, sum, na.rm = TRUE)
    temp <- sqrt((temp[2] + temp[3])/(temp[1] - 1))
    attr(temp, "names") <- NULL
    temp
}

RPSVx <- function(cx) {
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

autosigma <- function(capthist = NULL, bincapt, traps, mask, sv = NULL){
    obsRPSV <- RPSV.mod(bincapt, traps)
    secr:::naivesigma(obsRPSV, traps, mask, 0, 1)
}

autoD <- function(capthist = NULL, bincapt, traps, mask, sv){
    n <- dim(bincapt)[1]
    A <- attr(mask, "area")
    n/(sum(pdot(mask, traps, 0,
                     list(g0 = sv[2], sigma = sv[3]), 1))*A)
}

autog0 <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL){
    0.95
}

autosigmatoa <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL){
    0.0025
}

autokappa <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL){
    10
}

autossb0 <- function(capthist = NULL, bincapt, traps = NULL, mask = NULL, sv = NULL){
    159.9514
}

autossb1 <- function(capthist = NULL, bincapt, traps = NULL, mask = NULL, sv = NULL){
    0
}

autosigmass <- function(capthist = NULL, bincapt, traps = NULL, mask = NULL, sv = NULL){
    4.095078
}
