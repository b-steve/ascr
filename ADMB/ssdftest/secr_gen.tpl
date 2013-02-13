
DATA_SECTION

  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_matrix capt(1,345,1,6)
  init_matrix dist(1,6,1,3449)
  init_number trace

PARAMETER_SECTION

  objective_function_value f
  init_bounded_number D(0,1e+08)
  init_bounded_number par0(-50,50)
  init_bounded_number par1(-10,0)

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
  //p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;
  // Probability of detection at any trap for each location.
  for(i=1; i<=nmask; i++){
    p=1;
    for(j=1; j<=ntraps; j++){
      dvariable erfx = par0 - par1*dist(j,i);
      p1(j,i) = 0.5 - 0.5*(2*cumd_norm(erfx) - 1);
      p2(j,i) = 1 - p1(j,i);
      p*=p2(j)(i);
    }
    pm(i)=1-p;
  }
  logp1=log(p1 + DBL_MIN);
  logp2=log(p2 + DBL_MIN);
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
    cout << "D: " << D << ", par0: " << par0 << ", par1: " << par1 << ", loglik: " << -f << endl;
  }

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION
