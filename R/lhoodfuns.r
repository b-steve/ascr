## Likelihood function for simple SECR.
secrlikelihood <- function (beta, capthist, mask, dist = NULL, trace = FALSE) {

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

  L1 <- sum ( apply(capthist, 1, Dprwi,  prwi.s, log.gk, log.gk1, D) )


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

## Helper function for simple SECR.
Dprwi <- function (wi, prwi.s, log.gk, log.gk1, D) {
  ## wi is S x K, log.gk is K x M, prwi.s is S x M
  prwi.s <- wi %*% log.gk + (1-wi) %*% log.gk1
  ## sum log(p) over occasions, result is M-vector
  prwi.s <- apply(prwi.s,2,sum)
  log (sum(D * exp(prwi.s)))
}

## Likelihood function for SECR with TOA information.
secrlikelihood.toa1 <- function (beta, capthist, mask, dists=NULL, ssqtoa=NULL, trace=FALSE) {
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

## Helper functions for time of arrival analysis.
log.ftoa=function(capthist,ssqtoa,sigma.toa) {
  num.nonzeros=function(x) sum(x>0)
  logftoa=ssqtoa*0 # initialize
  M=apply(capthist,1,num.nonzeros)
  M2plus=which(M>1)
  madd=matrix(((1-M[M2plus])*log(sigma.toa)),nrow=dim(logftoa)[1],ncol=length(M2plus),byrow=TRUE)
  logftoa[,M2plus]=ssqtoa[,M2plus]/(-2*sigma.toa^2) + madd
  return(logftoa)
}

toa.ssq=function(wit,dists) {
  ssq=function(x) sum((x-mean(x))^2)
  v=330 # speed of sound
  wit.na=wit; wit.na[wit==0]=NA # mark those with no capture
  delt=na.omit(as.vector(wit.na)-dists/v) # need vector else get non-conformable arrays
  toassq=apply(delt,2,ssq)
  return(toassq)
}

## Likelihood function for SECR with angle information.
secrlikelihood.angs <- function (beta, capthist, mask, dists=NULL, angs=NULL, trace=FALSE) {
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

## Helper functions for SECR with angles.
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

## Likelihood function for SECR with signal strength information.
secrlikelihood.ss <- function (beta, capthist, mask, dists=NULL, cutoff, trace=FALSE) {
  ## 'real' parameter values
  D <- exp(beta[1])
  ## if D is modelled, expand here to a vector of length
  ## equal to number of rows in mask
  b0ss <- beta[2]
  b1ss <- beta[3]
  sigmass <- exp(beta[4]) # std. dev. of arrival time measurement error
  
  n <- nrow(capthist)        ## number observed
  S <- ncol(capthist)        ## number of occasions
  A <- attr( mask, 'area')   ## area of one cell
  
  ## distances if not passed as argument
  if (is.null(dists)) {
    traps <- traps(capthist)
    dists <- distances(traps, mask)
  }
  
  bincapt <- capthist
  bincapt[bincapt > 0] <- 1
  ## precompute probabilities for each detector and each mask point
  ## gk is K x M matrix
  muss <- b0ss + b1ss*dists
  gk <- 1 - pnorm(cutoff, muss, sigmass)
  log.gk <- log(gk + .Machine$double.xmin)
  log.gk1 <- log(1-gk + .Machine$double.xmin)
  
  ## probability of being caught at least once if at mask site m
  p.m <- 1 - apply(1-gk, 2, prod) ^ S  ## vector of length M
  sumDp <- sum(p.m * D)                ## scalar
  
  L1.X <- apply(capthist,1,log.Dprwi.ss, D, muss, sigmass, log.gk1)
  L1 <- sum(log(apply(exp(L1.X),2,sum)+.Machine$double.xmin))
  L2 <- - n * log(sumDp)
  L3 <- dpois (n, sumDp * A, log = TRUE)
  LL <- L1+L2+L3
  if (trace) {
    cat(D, b0ss, b1ss, sigmass, LL, '\n')
    flush.console()
  }
  if (!is.finite(LL))
    1e10
  else
    -LL   ## return negative log likelihood
}

## Helper function for SECR with angles.
log.Dprwi.ss <- function (wi, D, muss, sigmass, log.gk1) {
  ci <- wi
  ci[ci > 0] <- 1
  dens <- matrix(0, nrow = nrow(traps), ncol = nrow(mask))
  for (i in 1:nrow(traps)){
    dens[i, ] <- dnorm(wi[i], muss[i, ], sigmass, log = TRUE)
  }
  prwi.s <- ci %*% dens + (1-ci) %*% log.gk1
  prwi.s <-apply(prwi.s,2,sum)
  log(D) + prwi.s
}


## Darren's C++ implementation of SECR with angle information.
code <- '

	NumericVector BETA(beta);
	double D     = exp(BETA[0]);
        double G0    = 1/(1+exp(-BETA[1]));
	double SIGMA = exp(BETA[2]);

	int n    = as<int>(ncues);
	int K    = as<int>(ntraps);
	int M    = as<int>(npoints);

	int i,k,m;

	NumericMatrix RADS(radians);
	IntegerMatrix HASH1(hash1);
	IntegerMatrix HASH0(hash0);
	NumericMatrix DISTS(mask_dists);
	NumericMatrix ANGS(mask_angs);
	double A = as<double>(mask_area);
	
	NumericVector GK(K*M);
	NumericVector PM(M);
	PM.fill(1);
	for (m=0; m<M; m++) {
		for (k=0; k<K; k++) {
			GK(k*M+m) = G0*exp(-pow(DISTS(k,m),2) / 2 / pow(SIGMA,2));
			PM(m) *= 1-GK(k*M+m);
		}
	}
	double SUMDP = D*sum(1-PM);
	
	NumericMatrix L1_X(M,n);
	L1_X.fill(log(D));
	NumericVector LOG_GK = log(GK);
	int METHOD = as<int>(method);
	if(METHOD==1){
		double KAPPA = exp(BETA[3]);
		double CONSTANT = log(2*3.14159265359*Rf_bessel_i(KAPPA,0.0,1.0));
		for (int det=0; det<HASH1.nrow(); det++) {
			i = HASH1(det,0);
			k = HASH1(det,1);
			for (m=0; m<M; m++) {
				L1_X(m,i) += KAPPA*cos(RADS(i,k)-ANGS(k,m))-CONSTANT+LOG_GK(k*M+m);
			}
		}
	}else{
		for (int det=0; det<HASH1.nrow(); det++) {
			i = HASH1(det,0);
			k = HASH1(det,1);
			for (m=0; m<M; m++) {
				L1_X(m,i) += LOG_GK(k*M+m);
			}
		}
	}
	
	NumericVector LOG_GK1 = log(1-GK);
	for (int det=0; det<HASH0.nrow(); det++) {
		i = HASH0(det,0);
		k = HASH0(det,1);
		for (m=0; m<M; m++) {
			L1_X(m,i) += LOG_GK1(k*M+m);
		}
	}

	NumericVector L1(n);
	for (i=0; i<n; i++) {
		for (m=0; m<M; m++) {
			L1(i) += exp(L1_X(m,i));
		}
	}
	
	double NEG_LOGLIK = -(sum(log(L1))- n*log(SUMDP) + Rf_dpois(n,A*SUMDP,1));

	if(ISNAN(NEG_LOGLIK) | !finite(NEG_LOGLIK)){
		return wrap(1e10);
	}else{
		return wrap(NEG_LOGLIK);
	}
	'
	
secrlikelihood.cpp <- cxxfunction(signature(beta="numeric",
                                            method="integer",
                                            ncues="integer",
                                            ntraps="integer",
                                            npoints="integer",
                                            radians="numeric",
                                            hash1="integer",
                                            hash0="integer",
                                            mask_area="numeric",
                                            mask_dists="numeric",
                                            mask_angs="numeric") , body=code, plugin = "Rcpp")

