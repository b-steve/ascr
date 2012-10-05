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
  int i;
  f=0;
  for(i=1;i<=n;i++){
    g_cluster(i, D, g0, sigma, X(i), Y(i));
  }

SEPARABLE_FUNCTION void g_cluster(int i, const dvariable& D, const dvariable& g0, const dvariable& sigma, const dvariable& Xi, const dvariable& Yi)
  // Setting up variables
  int j,k;
  dvariable d,p,lambda,L1,L2,L3,indp,intL;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  if (i==1){
    // Probabilities of caputure at each location for each trap.
    // Add a small amount to prevent zeros.
    p1=g0*mfexp(-square(dist)*pow(2*square(sigma),-1))+DBL_MIN; // Why doesn't -square(dist)/(2*square(sigma)) work?!
    p2=1-p1;
    // Probability of detection at any trap for each location.
    // Used for numerical integration of *detectable* animal density.
    for(j=1; j<=nmask; j++){
      p=1;
      for(k=1; k<=ntraps; k++){
        p*=p2(k)(j);
      }
      pm(j)=1-p;
    }
    lambda=A*D*sum(pm);
    L2=-n*log(D*sum(pm));
    L3=log_density_poisson(n,lambda);
    f-=(L2+L3);
  }
  // Probability of capture histories for each animal.
  wi1=row(capt,i);
  // Prior on animal location.
  L1=log(D);
  for(j=1; j<=ntraps; j++){
    // Calculating distance between animal i and trap j.
    d=sqrt(square(Xi-trapsX(j))+square(Yi-trapsY(j)));
    // Calculating probability of capture of animal i at trap j.
    indp=g0*mfexp(-square(d)*pow(2*square(sigma),-1))+DBL_MIN;
    // Appropriate contribution to likelihood of animal capture/evasion.
    if (value(wi1(j))==1){
      L1+=log(indp); 
    }else{
      L1+=log(1-indp);
    }
  }
  f-=L1;

GLOBALS_SECTION
  #include <float.h>

