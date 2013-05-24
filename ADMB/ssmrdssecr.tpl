
PROCEDURE_SECTION
  // Setting up variables.
  const double pi=3.14159265359;
  int i, j;
  dvariable p, d, indivmuss;
  dvar_vector pm(1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix muss(1,ntraps,1,nmask);
  muss = //@LINKFN(ssb0 - ssb1*dist);
  for (i = 1; i <= nmask; i++){
    p = 1;
    for(j = 1; j <= ntraps; j++){
      p2(j,i) = cumd_norm((c - muss(j,i))/sigmass);
      p *= p2(j,i);
    }
    pm(i) = 1 - p + DBL_MIN;
  }
  dvariable L1 = 0;
  // Probability of capture histories for each animal.
  for (i = 1; i <= n; i++){
    for (j = 1; j <= ntraps; j++){
      d = indivdist(i,j);
      indivmuss = //@LINKFN(ssb0 - ssb1*d);
      if (capt(i,j) == 0){
        L1 += log(cumd_norm((c - indivmuss)/sigmass) + DBL_MIN);
      } else {
        L1 += -log(sigmass) - log(sqrt(2*pi)) + (square(sscapt(i,j) - indivmuss)/(-2*square(sigmass)));
      }
    }
  }
  L1 += n*log(D);
  dvariable L2 = -n*log(D*sum(pm));
  dvariable lambda = A*D*sum(pm);
  dvariable L3 = log_density_poisson(n,lambda);
  f = -(L1+L2+L3);
  if (trace == 1){
    //@TRACE;
  }

GLOBALS_SECTION
  #include <float.h>

