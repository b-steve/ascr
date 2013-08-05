
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
  random_effects_vector X(1,n)
  random_effects_vector Y(1,n)

PROCEDURE_SECTION
  int i;
  f=0;
  for(i=1;i<=n;i++){
    g_cluster(i, D, ssb0, ssb1, sigmass, X(i), Y(i));
  }

SEPARABLE_FUNCTION void g_cluster(int i, const dvariable& D, const dvariable& ssb0, const dvariable& ssb1, const dvariable& sigmass, const dvariable& Xi, const dvariable& Yi)
  // Setting up variables
  int j,k;
  const double pi=3.14159265359;
  dvariable d,p,lambda,L1,L2,L3,exprss;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix muss(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector ci1(1,ntraps);
  if (i==1){
    // Probabilities of caputure at each location for each trap.
    muss=ssb0+ssb1*dist;
    // Probability of detection at any trap for each location.
    // Used for numerical integration of *detectable* animal density.
    for(j=1; j<=nmask; j++){
      p=1;
      for(k=1; k<=ntraps; k++){
        p2(k,j)=cumd_norm((c-muss(k,j))/sigmass);
        p*=p2(k,j);
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
  ci1=row(sscapt,i);
  // Prior on animal location.
  L1=log(D);
  for(j=1; j<=ntraps; j++){
    // Calculating distance between animal i and trap j.
    d=sqrt(square(Xi-trapsX(j))+square(Yi-trapsY(j)));
    // Expected signal strength at trap.
    exprss=ssb0+ssb1*d;
    if (value(wi1(j))==1){
      // If detected, add log of normal density of observed received signal strength.
      L1+=-log(sigmass)-log(sqrt(2*pi))-(square(ci1(j)-exprss)/(2*square(sigmass)));
    }else{
      // If not detected, add log of probability of non-detection
      // (i.e. prob of signal strength being lower than the cutoff
      // (c)).
      L1+=log(cumd_norm((c-exprss)/sigmass)+DBL_MIN);
    }
  }
  f-=L1;

GLOBALS_SECTION
  #include <float.h>

TOP_OF_MAIN_SECTION
  arrmblsize=150000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_MAX_NVAR_OFFSET(4600404);
