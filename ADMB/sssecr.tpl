TOP_OF_MAIN_SECTION
  arrmblsize=10000000;

DATA_SECTION
  matrix capt(1,n,1,ntraps)
  !!for(int i=1;i<=n;i++){
  !!  for(int j=1;j<=ntraps;j++){
  !!    if(sscapt(i,j)>0){
  !!      capt(i,j)=1;
  !!    }
  !!    else{
  !!      capt(i,j)=0;
  !!    }
  !!  }
  !!}

PROCEDURE_SECTION
  // Setting up variables
  const double pi=3.14159265359;
  int i,j;
  dvariable p,lambda,L1,L2,L3;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);
  dvar_vector ssll(1,nmask);
  dvar_vector ess(1,nmask);
  // Probabilities of caputure at each location for each trap.
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
  }
  L1=0;
  // Probability of capture histories for each animal.
  for(i=1; i<=n; i++){
    wi1=capt(i)(1,ntraps);
    wi2=1-wi1;
    ssll=0;
    // Likelihood due to signal strengths.
    for(j=1; j<=ntraps; j++){
      if (capt(i)(j)==1){
        ess=ssb0+ssb1*row(dist,j);
        ssll+=-log(sigmass)-(square(sscapt(i)(j)-ess)/(2*square(sigmass)));
      }
    }
    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+ssll))+DBL_MIN);
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION
