
PROCEDURE_SECTION
  int i,j;
  dvariable p,lambda,L1,L2,L3;
  dvar_matrix maskp1(1,ntraps,1,nmask);
  dvar_matrix maskp2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_matrix indivp1(1,n,1,ntraps);
  dvar_matrix indivp2(1,n,1,ntraps);
  dvar_matrix logindivp1(1,n,1,ntraps);
  dvar_matrix logindivp2(1,n,1,ntraps);
  dvar_matrix logprobs(1,n,1,ntraps);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);
  // Probability of capture for each individual at each trap.
  maskp1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;
  maskp2=1-maskp1;
  logp1=log(maskp1);
  logp2=log(maskp2);
  // Probability of detection at any trap for each location.
  for(i=1; i<=nmask; i++){
    p=1;
    for(j=1; j<=ntraps; j++){
      p*=maskp2(j)(i);
    }
    pm(i)=1-p;
  }
  // Probability of capture histories for each animal.
  indivp1=g0*mfexp(-square(indivdist)/(2*square(sigma)))+DBL_MIN;
  indivp2=1-indivp1;
  logindivp1=log(indivp1);
  logindivp2=log(indivp2);
  logprobs=elem_prod(logindivp1,capt)+elem_prod(logindivp2,1-capt);
  L1=sum(logprobs)+n*log(D);
  L2=-n*log(D*sum(pm));
  lambda=A*D*sum(pm);
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  if (trace == 1){
  cout << "D: " << D << ", g0: " << g0 << ", sigma: " << sigma << ", loglik: " << -f << endl;
  }

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION

