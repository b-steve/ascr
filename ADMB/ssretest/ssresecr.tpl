TOP_OF_MAIN_SECTION
  arrmblsize=2500000;

DATA_SECTION
  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_number c
  //init_matrix sscapt(1,345,1,6)
  //init_matrix dist(1,6,1,3449)
  //init_matrix capt(1,345,1,6)
  init_matrix sscapt(1,25,1,6)
  init_matrix dist(1,6,1,3449)
  init_matrix capt(1,25,1,6)
  init_matrix traps(1,6,1,2)
  init_number trace
  vector trapsX(1,ntraps)
  vector trapsY(1,ntraps)
  !!trapsX=column(traps,1);
  !!trapsY=column(traps,2);

PARAMETER_SECTION

  objective_function_value f
  init_bounded_number D(0,1e+07)
  init_number ssb0mu
  init_bounded_number ssb0sigma(0,1e+05)
  init_bounded_number ssb1(-10,-0.001)
  init_bounded_number sigmass(0,1e+05)
  random_effects_vector X(1,n)
  random_effects_vector Y(1,n)
  random_effects_vector ssb0(1,n)

PROCEDURE_SECTION
  int i;
  f=0;
  for(i=1;i<=n;i++){
    g_cluster(i, D, ssb0mu, ssb0sigma, ssb1, sigmass, X(i), Y(i), ssb0(i));
  }
  cout << "D: " << D << " ssb0mu: " << ssb0mu << " ssb0sigma: " << ssb0sigma << " ssb1: " << ssb1 << " sigmass: " << sigmass << endl;

SEPARABLE_FUNCTION void g_cluster(int i, const dvariable& D, const dvariable& ssb0mu, const dvariable& ssb0sigma, const dvariable& ssb1, const dvariable& sigmass, const dvariable& Xi, const dvariable& Yi, const dvariable& ssb0i)
  // Setting up variables
  int j,k;
  const double pi=3.14159265359;
  dvariable d,p,lambda,L1,L2,L3,indp,expss,sigmarec;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix muss(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector ci1(1,ntraps);
  if (i==1){
    muss=ssb0mu+ssb1*dist;
    sigmarec=sigmass;
    sigmarec=pow(square(ssb0sigma)+square(sigmass),0.5);
    // Probability of detection at any trap for each location.
    // Used for numerical integration of *detectable* animal density.
    for(j=1; j<=nmask; j++){
      p=1;
      for(k=1; k<=ntraps; k++){
        p2(k,j)=cumd_norm((c-muss(k,j))/sigmarec);
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
  ci1=row(capt,i);
  wi1=row(sscapt,i);
  // Prior on animal location.
  L1=log(D);
  // Prior on animal signal strength.
  L1+=-log(ssb0sigma)-log(sqrt(2*pi))+(square(ssb0i)/(-2*square(ssb0sigma)));
  for(j=1; j<=ntraps; j++){
    // Calculating distance between animal i and trap j.
    d=sqrt(square(Xi-trapsX(j))+square(Yi-trapsY(j)));
    // Calculating probability of capture of animal i at trap j.
    expss=ssb0mu+ssb0i+ssb1*d;
    // Appropriate contribution to likelihood of animal capture/evasion.
    if (value(ci1(j))==1){
      L1+=-log(sigmass)-log(sqrt(2*pi))+(square(wi1(j)-expss)/(-2*square(sigmass)));
    }else{
      L1+=log(cumd_norm((c-expss)/sigmass)+DBL_MIN);
    }
  }
  f-=L1;

GLOBALS_SECTION
  #include <float.h>

