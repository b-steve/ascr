
PROCEDURE_SECTION
  int i;
  dvariable p,lambda,L1,L2,L3;
  dvar_vector maskp11(1,nmask);
  dvar_vector maskp12(1,nmask);
  dvar_vector maskp21(1,nmask);
  dvar_vector maskp22(1,nmask);
  dvar_vector indivprobs(1,n);
  dvar_vector logprobs(1,n);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);
  // Probability of capture for each individual at each trap.
  for (i=1; i<=nmask; i++){
    maskp11(i)=g01*mfexp(-square(dist(1,i))/(2*square(sigma1)))+DBL_MIN;
    maskp21(i)=g02*mfexp(-square(dist(2,i))/(2*square(sigma2)))+DBL_MIN;
  }
  maskp12=1-maskp11;
  maskp22=1-maskp21;
  // Probability of detection at any trap for each location.
  for(i=1; i<=nmask; i++){
    p=1;
    p*=maskp12(i);
    p*=maskp22(i);
    pm(i)=1-p;
  }
  // Probability of capture histories for each animal.
  for (i=1; i<=n; i++){
    indivprobs(i)=g01*mfexp(-square(indivdist(i,1))/(2*square(sigma1)))+DBL_MIN;
  }
  logprobs=elem_prod(log(indivprobs),column(capt,1))+elem_prod(log(1-indivprobs),column(1-capt,1));
  for (i=1; i<=n; i++){
    indivprobs(i)=g02*mfexp(-square(indivdist(i,1))/(2*square(sigma2)))+DBL_MIN;
  }
  logprobs+=elem_prod(log(indivprobs),column(capt,2))+elem_prod(log(1-indivprobs),column(1-capt,2));
  L1=sum(logprobs)+n*log(D);
  L2=-n*log(D*sum(pm));
  lambda=A*D*sum(pm);
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  if (trace == 1){
    cout << "D: " << D << ", g01: " << g01 << ", sigma1: " << sigma1 << ", g02: " << g02 << ", sigma2: " << sigma2 << ", loglik: " << -f << endl;
  }

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION

