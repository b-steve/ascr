  sdreport_number esa
// Flag for creating sdreport_number for D
//@SDREPD

PROCEDURE_SECTION
  // Setting up variables
  const double DBL_MIN = 1e-150;
  const double pi = 3.14159265359;
  int i, j;
  dvariable p, d;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  f = 0;
  // Flag for specifying D
  //@SPECD
  // Probability of detection at any trap for each location.
  for (i = 1; i <= nmask; i++){
    p = 1;
    for (j = 1; j <= ntraps; j++){
      d = dist(j,i);
      // Flag for detection function insertion.
      p1(j,i) = //@DETFN;
      p2(j,i) = 1 - p1(j,i);
      p *= p2(j,i);
    }
    pm(i) = 1 - p + DBL_MIN;
  }
  logp1 = log(p1 + DBL_MIN);
  logp2 = log(p2 + DBL_MIN);
  dvariable L1 = 0;
  // Probability of capture histories for each animal.
  for (i = 1; i <= n; i++){
    wi1 = row(capt,i);
    L1 += log(D*sum(mfexp(wi1*logp1 + (1 - wi1)*logp2)) + DBL_MIN);
  }
  // Calculating esa.
  esa = A*sum(pm);
  // Calculating esa.
  esa = A*sum(pm);
  // Putting log-likelihood together.
  dvariable lambda = D*esa;
  dvariable L2 = -n*log(D*sum(pm));
  dvariable L3 = log_density_poisson(n,lambda);
  f -= L1 + L2 + L3;
  if (trace == 1){
    //@TRACE;
  }



