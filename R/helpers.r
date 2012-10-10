## Miscellaneous helpful functions.

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
distances.cpp <- cxxfunction(signature(traps = "numeric", mask = "numeric"),
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
angles.cpp <- cxxfunction(signature(traps = "numeric", mask = "numeric") ,
                           body = angcode, plugin = "Rcpp")

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

## David's function to determine frog ID of a click.
make.frog.captures=function(mics,clicks,dt){
  K=dim(mics)[1]
  captures=clicks
  ## Times of clicks in current set.
  ct=rep(-Inf,K)
  ## Counter for click.
  ID=1
  ## Store time of click by mic 1 in ct[mic].
  ct[clicks$trap[1]]=clicks$tim[1]
  ## Indicator that is true when click can't be part of current set.
  new=FALSE
  nclicks=length(clicks$tim)
  for(i in 2:nclicks){
    if(ct[clicks$trap[i]]>-Inf){
      nd=length(which(ct>-Inf))
      ## Make all but last of those in current set part of same capture history
      captures$ID[(i-nd):(i-1)]=ID
      ## Re-initialise
      ct=rep(-Inf,K)
      ## Store time of click by mic i in ct[mic]
      ct[clicks$trap[i]]=clicks$tim[i] 
      ID=ID+1
      if(i==nclicks) captures$ID[i]=ID # write last record with new ID      
    }
    else { # next click on a mic not in current set
      ct[clicks$trap[i]]=clicks$tim[i] # store time of click by mic i in ct[mic]
      ctset=which(ct>-Inf)
      dts=dt[ctset,clicks$trap[i]] # times between mics in current set and mic i
      cts=-(ct[ctset]-clicks$tim[i]) # times between clicks in current set and click i
      if(any((cts-dts)>0)) new=TRUE
      if(new) { 
        nd=length(which(ct>-Inf))-1 # number clicks in current set before new click added
        captures$ID[(i-nd):(i-1)]=ID # make all but last of those in current set part of same capture history
        ct=rep(-Inf,K) # re-initialise 
        ct[clicks$trap[i]]=clicks$tim[i] # store time of click by mic i in ct[mic]
        ID=ID+1
        new=FALSE
        if(i==nclicks) captures$ID[i]=ID # write last record with new ID
      } else if(i==nclicks){
        nd=length(which(ct>-Inf)) # number clicks in current set
        captures$ID[(i-nd+1):i]=ID # make all in current set part of same capture history      
      }
    }
  }
  return(captures)
}


