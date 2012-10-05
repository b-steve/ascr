
TOP_OF_MAIN_SECTION
  arrmblsize=20000000;

PROCEDURE_SECTION
  // Setting up variables
  int i,j;
  dvariable p,lambda,L1,L2,L3;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);
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
    wi1=row(capt,i);
    wi2=1-wi1;
    L1+=log(D*sum(mfexp(wi1*logp1+wi2*logp2)));
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  if (trace == 1){
    cout << "D: " << D << ", g0: " << g0 << ", sigma: " << sigma << ", loglik: " << -f << endl;
  }

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION
