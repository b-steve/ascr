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

autoini.mod <- function(capthist, trps, mask, detectfn = 0, thin = 0.2) {
  naivedcall <- function(sigma) {
    temp <- .C("naived", PACKAGE = "secr", as.double(sigma), 
               as.integer(k), as.integer(m), as.double(unlist(trps)), 
               as.double(unlist(mask)), as.integer(detectfn), value = double(1))
    db - temp$value
  }
  naivecap2 <- function(g0, sigma, cap) {
    temp <- .C("naivecap2", PACKAGE = "secr", as.integer(prox), 
               as.double(g0), as.double(sigma), as.integer(s), as.integer(k), 
               as.integer(m), as.double(unlist(trps)), as.double(unlist(mask)), 
               as.integer(detectfn), value = double(1))
    cap - temp$value
  }
  naiveesa <- function(g0, sigma) {
    nc <- 1
    g0sigma0 <- matrix(rep(c(g0, sigma), c(2, 2)), nrow = 2)
    gs0 <- rep(1, 2 * s * k)
    area <- attr(mask, "area")
    param <- 0
    miscparm <- 1
    temp <- try(.C("integralprw1", PACKAGE = "secr", as.integer(dettype), 
                   as.integer(param), as.double(g0sigma0), as.integer(nc), 
                   as.integer(s), as.integer(k), as.integer(m), as.integer(1), 
                   as.double(unlist(trps)), as.double(unlist(mask)), 
                   as.integer(nrow(g0sigma0)), as.integer(gs0), as.integer(1), 
                   as.double(area), as.double(miscparm), as.integer(detectfn), 
                   as.integer(0), as.integer(0), a = double(nc), resultcode = integer(1)))
    if (temp$resultcode != 0) 
      stop("error in external function 'integralprw1'; ", 
           "possibly the mask is too large")
    temp$a
  }
  if (nrow(capthist) < 5) 
    stop("too few values in session 1 to determine start; set manually")
  if (is.character(detectfn)) 
    detectfn <- detectionfunctionnumber(detectfn)
  if (!(detectfn %in% c(0))) 
    stop("only halfnormal detection function implemented in 'autoini'")
  prox <- detector(trps) %in% c("proximity", "count", "signal", 
                                "times")
  n <- nrow(capthist)
  s <- ncol(capthist)
  k <- nrow(trps)
  if ((nrow(mask) > 100) & (thin > 0) & (thin < 1)) 
    mask <- mask[runif(nrow(mask)) < thin, ]
  else thin <- 1
  m <- nrow(mask)
  if (length(dim(capthist)) > 2) 
    cpa <- sum(abs(capthist))/n
  else cpa <- sum(abs(capthist) > 0)/n
  obsRPSV <- RPSV(capthist)
  if (is.na(obsRPSV) | (obsRPSV < 1e-10)) {
    db <- dbar(capthist)
    if (!is.null(attr(trps, "spacing"))) {
      if (is.na(db)) {
        warning("could not calculate 'dbar'; using detector spacing")
        db <- attr(trps, "spacing")
      }
      if (db < (attr(trps, "spacing")/100)) {
        warning("'dbar' close to zero; using detector spacing instead")
        db <- attr(trps, "spacing")
      }
    }
    if (is.na(db) | is.nan(db) | (db < 1e+10)) 
      return(list(D = NA, g0 = NA, sigma = NA))
    else tempsigma <- uniroot(naivedcall, lower = db/10, 
                              upper = db * 10, tol = 0.001)$root
  }
  else {
    tempsigma <- naivesigma(obsRPSV = obsRPSV, trps = trps, 
                            mask = mask, detectfn = detectfn, z = 1)
  }
  low <- naivecap2(0.001, sigma = tempsigma, cap = cpa)
  upp <- naivecap2(0.999, sigma = tempsigma, cap = cpa)
  badinput <- FALSE
  if (is.na(low) | is.na(upp)) 
    badinput <- TRUE
  else if (sign(low) == sign(upp)) 
    badinput <- TRUE
  if (badinput) {
    warning("'autoini' failed to find g0; setting initial g0 = 0.1")
    tempg0 <- 0.1
  }
  else tempg0 <- uniroot(naivecap2, lower = 0.001, upper = 0.999, 
                         f.lower = low, f.upper = upp, tol = 0.001, sigma = tempsigma, 
                         cap = cpa)$root
  esa <- naiveesa(tempg0, tempsigma)
  list(D = n/esa * thin, g0 = tempg0, sigma = tempsigma)
}
