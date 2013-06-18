
// Flag for creating sdreport_number for D
//@SDREPD

PROCEDURE_SECTION
  // Setting up variables.
  const double DBL_MIN = 1e-150;
  int i,j;
  dvariable p, p1, d;
  dvar_matrix indivp1(1,n,1,ntraps);
  dvar_matrix indivp2(1,n,1,ntraps);
  dvar_matrix logindivp1(1,n,1,ntraps);
  dvar_matrix logindivp2(1,n,1,ntraps);
  dvar_matrix logprobs(1,n,1,ntraps);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);
  f = 0;
  // Flag for specifying D
  //@SPECD
  // Probability of detection at any trap for each location.
  for (i = 1; i <= nmask; i++){
    p = 1;
    for(j = 1; j <= ntraps; j++){
      d = dist(j,i);
      // Flag for detection function insertion.
      p1 = //@DETFN;
      p *= 1 - p1;
    }
    pm(i) = 1 - p + DBL_MIN;
  }
  // Probability of capture histories for each animal.
  for (i = 1; i <= n; i++){
    for (j = 1; j <= ntraps; j++){
      d = indivdist(i,j);
      indivp1(i,j) = //@DETFN;
      indivp2(i,j) = 1 - indivp1(i,j);
    }
  }
  logindivp1 = log(indivp1 + DBL_MIN);
  logindivp2 = log(indivp2 + DBL_MIN);
  logprobs = elem_prod(logindivp1,capt) + elem_prod(logindivp2,1 - capt);
  dvariable L1 = sum(logprobs) + n*log(D);
  dvariable L2 = -n*log(D*sum(pm));
  dvariable lambda = A*D*sum(pm);
  dvariable L3 = log_density_poisson(n,lambda);
  f -= L1 + L2 + L3;
  if (trace == 1){
    //@TRACE;
  }


