DATA_SECTION
  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_matrix capt(1,40,1,6)
  init_matrix sscapt(1,40,1,6)
  init_matrix dist(1,6,1,3449)
  init_number c
  init_matrix traps(1,6,1,2)
  init_number trace
  vector trapsX(1,ntraps)
  vector trapsY(1,ntraps)
  !!trapsX=column(traps,1);
  !!trapsY=column(traps,2);

PARAMETER_SECTION
  objective_function_value f
  init_bounded_number D(0,1e+07)
  init_number ssb0
  init_bounded_number ssb1(-10,0)
  init_bounded_number sigmass(0,1e+07)

PROCEDURE_SECTION
  // Setting up variables
  const double pi=3.14159265359;
  int i,j,k;
  dvariable p,lambda,L1,L2,L3;
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector ci1(1,ntraps);
  dvar_matrix muss(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  //muss=mfexp(ssb0+ssb1*dist);
  muss=ssb0+ssb1*dist;
  for(i=1; i<=nmask; i++){
    p=1;
    for(j=1; j<=ntraps; j++){
      p2(j,i)=cumd_norm((c-muss(j,i))/sigmass);
      p*=p2(j,i);
    }
    pm(i)=1-p;
  }
  // Probability of non-detection.
  logp2=log(p2+DBL_MIN);
  L1=0;
  // Probability of capture histories for each animal.
  for(i=1; i<=n; i++){
    logp1=0;
    wi1=row(sscapt,i);
    ci1=row(capt,i);
    for(j=1; j<=ntraps; j++){
      if (ci1(j)==1){
        // If animal detected, calculate log of normal density of observed received signal strength.
        logp1(j)(1,nmask)=-log(sigmass)-log(sqrt(2*pi))+(square(wi1(j)-row(muss,j))/(-2*square(sigmass)));
      }
    }
    L1+=log(D*sum(mfexp(ci1*logp1+(1-ci1)*logp2)+DBL_MIN));
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm)+DBL_MIN;
  L2=-n*log(D*sum(pm)+DBL_MIN);
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  if (trace==1){
    cout << "D: " << D << ", ssb0: " << ssb0 << ", ssb1: " << ssb1 << ", sigmass: " << sigmass << ", loglik: " << -f << endl;
  }


GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION
