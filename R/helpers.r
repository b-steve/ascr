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

distcode <- '
  NumericMatrix TRAPS(traps);
  NumericMatrix MASK(mask);
  int K = TRAPS.nrow();
  int M = MASK.nrow();
  NumericMatrix DISTANCES(K,M);
  for (int m=0; m<M; m++) {
  	for (int k=0; k<K; k++) {
      DISTANCES(k,m) = pow(pow(TRAPS(k,0)-MASK(m,0),2)+pow(TRAPS(k,1)-MASK(m,1),2),0.5);
    }
  }
  return wrap(DISTANCES);
'
distances.cpp <-  cxxfunction(signature(traps = "numeric", mask = "numeric"),
                              body = distcode, plugin = "Rcpp")

angles <- function (X, Y) {
# X and Y are 2-column matrices of coordinates Returns angles (0,360]
# from points in X to points in Y in matrix of dimension nrow(X) x
# nrow(Y)
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

angcode <- '
  NumericMatrix TRAPS(traps);
	NumericMatrix MASK(mask);
	int K = TRAPS.nrow();
	int M = MASK.nrow();
	NumericMatrix ANGLES(K,M);
	double X;
	double Y;
	double pi = 3.14159265359;
	for (int k=0; k<K; k++) {
		for (int m=0; m<M; m++) {
			X = MASK(m,0)-TRAPS(k,0);
			Y = MASK(m,1)-TRAPS(k,1);
			ANGLES(k,m) = atan(X/Y);
			if(Y<0) ANGLES(k,m) += pi;
			if(X<0 & Y>=0) ANGLES(k,m) += 2*pi;
		}
	}
	return wrap(ANGLES);
	'
angles.cpp <-  cxxfunction(signature(traps = "numeric", mask = "numeric") ,
                           body = angcode, plugin = "Rcpp")


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

Dprwi <- function (wi, prwi.s, log.gk, log.gk1, D) {
  ## wi is S x K, log.gk is K x M, prwi.s is S x M
  prwi.s <- wi %*% log.gk + (1-wi) %*% log.gk1
  ## sum log(p) over occasions, result is M-vector
  prwi.s <- apply(prwi.s,2,sum)
  log (sum(D * exp(prwi.s)))
}

## returns sum over detectors on which animal was detected, of squared
## (t0hat-mean(t0hat)), where t0hat_m=toa_m-dists_m/v, # toa_m is time
## of arrival at detector m, dists_m is distance from detector m to
## each mask point.
##
## Returns vecor of length=number of mask points
toa.ssq=function(wit,dists) {

  ssq=function(x) sum((x-mean(x))^2)
  v=330 # speed of sound
  wit.na=wit; wit.na[wit==0]=NA # mark those with no capture
  delt=na.omit(as.vector(wit.na)-dists/v) # need vector else get non-conformable arrays
  toassq=apply(delt,2,ssq)
  return(toassq)
}

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

log.ftoa=function(capthist,ssqtoa,sigma.toa) {
  num.nonzeros=function(x) sum(x>0)
  logftoa=ssqtoa*0 # initialize
  M=apply(capthist,1,num.nonzeros)
  M2plus=which(M>1)
  madd=matrix(((1-M[M2plus])*log(sigma.toa)),nrow=dim(logftoa)[1],ncol=length(M2plus),byrow=TRUE)
  logftoa[,M2plus]=ssqtoa[,M2plus]/(-2*sigma.toa^2) + madd
  return(logftoa)
}

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

## Following are helper functions for naive sigma estimation.
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

## The following are functions for automatic start value generation.
autosigma <- function(capthist = NULL, bincapt, traps, mask, sv = NULL, cutoff = NULL,
                      method = NULL){
  obsRPSV <- RPSV.mod(bincapt, traps)
  secr:::naivesigma(obsRPSV, traps, mask, 0, 1)
}

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

autog0 <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                   cutoff = NULL, method = NULL){
  0.95
}

autosigmatoa <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                         cutoff = NULL, method = NULL){
  0.0025
}

autokappa <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                      cutoff = NULL, method = NULL){
  10
}

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

autossb1 <- function(capthist = NULL, bincapt = NULL, traps = NULL, mask = NULL, sv = NULL,
                     cutoff = NULL, method = NULL){
  -0.1
}

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

## Darren's C++ stuff:
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

## Functions for making a .tpl.

make.top.of.main <- function(memory){
  if (!is.null(memory)){
    cat("TOP_OF_MAIN_SECTION\n  arrmblsize=", memory, ";", file = "secr.tpl", sep = "") 
  }
}

make.common.variable.init <- function(){
  string <- "\n\nPROCEDURE_SECTION
  // Setting up variables
  int i,j;
  dvariable p,lambda,L1,L2,L3;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);" 
  cat(string, file = "secr.tpl", sep = "", append = TRUE)
}

make.extra.variable.init <- function(methods){
  toastring <- "\n  dvariable nzz;
  dvar_vector toall(1,nmask);"["toa" %in% methods]
  angstring <- "\n  const double pi=3.14159265359;
  dvar_vector angll(1,nmask);"["ang" %in% methods]
  ssstring <- "\n  const double pi=3.14159265359;
  dvar_vector ssll(1,nmask);
  dvar_vector ess(1,nmask);"["ss" %in% methods]
  cat(toastring, angstring, ssstring, file = "secr.tpl", sep = "", append = TRUE)
}

make.probabilities <- function(){
  string <- "\n  // Probabilities of caputure at each location for each trap.
  // Add a small amount to prevent zeros.
  p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;
  p2=1-p1;
  logp1=log(p1);
  logp2=log(p2);
  // Probability of detection at any trap for each location.
  for(i=1; i<=nmask; i++){
    p=1;
    for(j=1; j<=ntraps; j++){
      p*=p2(j)(i);
    }
    pm(i)=1-p;
  }"
  cat(string, file = "secr.tpl", append = TRUE)
}

make.L1.likelihood <- function(methods){
  common.start <- "\n  L1=0;
  // Probability of capture histories for each animal.
  for(i=1; i<=n; i++){
    wi1=capt(i)(1,ntraps);
    wi2=1-wi1;"
  toastring <-
    "\n    nzz=sum(wi1);
    toall=(1-nzz)*log(sigmatoa)-((row(toassq, i))/(2*square(sigmatoa)));"["toa" %in% methods]
  angstring <-
    "\n    angll=0;
    // Likelihood due to angles.
    for(j=1; j<=ntraps; j++){
      // Von-Mises density contribution for each trap.
      if(capt(i)(j)==1){
        angll+=kappa*cos(angcapt(i)(j)-row(ang,j));
      }
    }
    // Term in Von-Mises density not dependent on data.
    angll-=sum(wi1)*log(2*pi*bessi0(kappa));"["ang" %in% methods]
  ssstring <-
    "\n    ssll=0;
    // Likelihood due to signal strengths.
    for(j=1; j<=ntraps; j++){
      if (capt(i)(j)==1){
        ess=ssb0+ssb1*row(dist,j);
        ssll+=-log(sigmass)-(square(sscapt(i)(j)-ess)/(2*square(sigmass)));
      }
    }"["ss" %in% methods]
  common.end <- paste("\n    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)",
                      "+toall"["toa" %in% methods], "+angll"["ang" %in% methods],
                      "+ssll"["ss" %in% methods], "))+DBL_MIN);\n  }", sep = "")
  cat(common.start, toastring, angstring, ssstring, common.end, file = "secr.tpl",
      sep = "", append = TRUE)
}

make.together.likelihood <- function(){
  string <- "\n  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);"
  cat(string, file = "secr.tpl", append = TRUE)
}

make.globals <- function(methods){
  common.string <- "\n\nGLOBALS_SECTION
  #include <float.h>"
  ss.string <- "\n  #include <bessel.cxx>"["ang" %in% methods]
  cat(common.string, ss.string, "\n\n", file = "secr.tpl", sep = "", append = TRUE)
}

make.all.tpl <- function(memory, methods){
  if (file.exists("secr.tpl")){
    file.remove("secr.tpl")
  }
  file.create("secr.tpl")
  make.top.of.main(memory)
  make.common.variable.init()
  make.extra.variable.init(methods)
  make.probabilities()
  make.L1.likelihood(methods)
  make.together.likelihood()
  make.globals(methods)
}

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
