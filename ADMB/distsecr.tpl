
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
  dvar_vector distll(1,nmask);
  dvar_vector beta(1,nmask);
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
    distll=0;
    // Likelihood due to distances.
    for(j=1; j<=ntraps; j++){
      // Gamma density contribution for each trap.
      if(capt(i)(j)==1){
	beta=alpha/row(dist,j);
	distll+=alpha*log(beta)+(alpha-1)*log(distcapt(i)(j))-(beta*distcapt(i)(j))-gammln(alpha);
      }
    }
    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+distll))+DBL_MIN);
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  if (trace == 1){
  cout << "D: " << D << ", g0: " << g0 << ", sigma: " << sigma << ", alpha: " << alpha << ", loglik: " << -f << endl;
  }

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION
