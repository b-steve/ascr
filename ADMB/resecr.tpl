TOP_OF_MAIN_SECTION
  arrmblsize=150000000;

DATA_SECTION

  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_matrix capt(1,20,1,6)
  init_matrix dist(1,6,1,3449)
  init_matrix traps(1,6,1,2)
  init_number trace
  vector trapsX(1,ntraps)
  vector trapsY(1,ntraps)
  !!trapsX=column(traps,1);
  !!trapsY=column(traps,2);

PARAMETER_SECTION

  objective_function_value f
  init_bounded_number D(0,1e+07)
  init_bounded_number g0(0,1)
  init_bounded_number sigma(0,1e+05)
  random_effects_vector X(1,n)
  random_effects_vector Y(1,n)

PROCEDURE_SECTION
  // Setting up variables
  int i,j;
  dvariable d,p,lambda,L1,L2,L3,indp,intL;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  // Probabilities of caputure at each location for each trap.
  // Add a small amount to prevent zeros.
  p1=g0*mfexp(-square(dist)*pow(2*square(sigma),-1))+DBL_MIN; // Why doesn't -square(dist)/(2*square(sigma)) work?!
  p2=1-p1;
  // Probability of detection at any trap for each location.
  // Used for numerical integration of *detectable* animal density.
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
    wi1=row(capt,i);
    // Prior on animal location.
    intL=log(D);
    for(j=1; j<=ntraps; j++){
      // Calculating distance between animal i and trap j.
      d=sqrt(square(X(i)-trapsX(j))+square(Y(i)-trapsY(j)));
      // Calculating probability of capture of animal i at trap j.
      indp=g0*mfexp(-square(d)*pow(2*square(sigma),-1))+DBL_MIN;
      // Appropriate contribution to likelihood of animal capture/evasion.
      if (value(wi1(j))==1){
        intL+=log(indp); 
      }else{
        intL+=log(1-indp);
      }
    }
    L1+=intL;
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);

GLOBALS_SECTION
  #include <float.h>

