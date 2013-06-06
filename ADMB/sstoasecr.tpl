
PROCEDURE_SECTION
  // Setting up variables
  const double pi=3.14159265359;
  int i, j;
  dvariable p, d, nzz;
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector ci1(1,ntraps);
  dvar_vector toall(1,nmask);
  dvar_matrix muss(1,ntraps,1,nmask);
  muss = //@LINKFN(ssb0 - ssb1*dist);
  for (i = 1; i <= nmask; i++){
    p=1;
    for (j = 1; j <= ntraps; j++){
      p2(j,i) = cumd_norm((c - muss(j,i))/sigmass);
      p *= p2(j,i);
    }
    pm(i) = 1 - p + DBL_MIN;
  }
  logp2 = log(p2 + DBL_MIN);
  dvariable L1=0;
  // Probability of capture histories for each animal.
  for (i = 1; i <= n; i++){
    logp1 = 0;
    wi1 = row(sscapt,i);
    ci1 = row(capt,i);
    nzz = sum(ci1);
    for (j = 1; j <= ntraps; j++){
      if (ci1(j) == 1){
        logp1(j)(1,nmask)= -log(sigmass) - log(sqrt(2*pi)) + (square(wi1(j)-row(muss,j))/(-2*square(sigmass)));
      }
    }
    toall = (1 - nzz)*log(sigmatoa) - ((row(toassq,i))/(2*square(sigmatoa)));
    L1 += log(D*sum(mfexp((ci1*logp1 + (1 - ci1)*logp2) + toall) + DBL_MIN));
  }
  // Putting log-likelihood together.
  dvariable lambda = A*D*sum(pm);
  dvariable L2 = -n*log(D*sum(pm));
  dvariable L3 = log_density_poisson(n,lambda);
  f = -(L1+L2+L3);
  if (trace==1){
    //@TRACE;
  }

GLOBALS_SECTION
  #include <float.h>


