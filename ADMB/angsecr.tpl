TOP_OF_MAIN_SECTION
  arrmblsize=10000000;

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
  dvar_vector angll(1,nmask);
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
    angll=0;
    // Likelihood due to angles.
    for(j=1; j<=ntraps; j++){
      // Von-Mises density contribution for each trap.
      if(capt(i)(j)==1){
        angll+=kappa*cos(angcapt(i)(j)-row(ang,j));
      }
    }
    // Term in Von-Mises density not dependent on data.
    angll-=sum(wi1)*log(2*pi*bessi0(kappa));
    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+angll)));
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);

GLOBALS_SECTION
  #include <float.h>
  #include <bessel.cxx>

REPORT_SECTION
