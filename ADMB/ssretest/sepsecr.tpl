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
  random_effects_vector m(1,2)

PROCEDURE_SECTION
  int i;
  f=0;
  ll_exp(D, g0, sigma, m(1), m(2));
  for(i=1;i<=n;i++){
    g_cluster(i, D, g0, sigma, X(i), Y(i));
  }

  // Exponential part of the likelihood.
SEPARABLE_FUNCTION void ll_exp(const dvariable& D, const dvariable& g0, const dvariable& sigma, const dvariable& m1, const dvariable& m2)
  const double pi = 3.14159265359;
  dvariable d, muss, p1;
  f -= -0.5*square(m1) - 0.5*log(2*pi);
  f -= -0.5*square(m2) - 0.5*log(2*pi);
  // Variables with uniform distribution on [-40, 40]
  dvariable x = 80.0*cumd_norm(m1) - 40.0;
  dvariable y = 80.0*cumd_norm(m2) - 40.0;
  dvariable p = 1.0;
  for (int j = 1; j <= ntraps; j++){
    d = sqrt(square(x - trapsX(j)) + square(y - trapsY(j)) + DBL_MIN);
    p1 = g0*mfexp(-square(d)*pow(2*square(sigma),-1)) + DBL_MIN;
    p *= 1 - p1;
  }
  // Probability of capture.
  dvariable pm = 1 - p;

  // Contribution to the likelihood is: exp(-1 * \int D * Area * Pr(capture) * f(X) dX)
  // Here we need to give: log(D * Area * Pr(capture) * f(X))
  f -= log(D);

  // Area of support of X = 80^2/10000 (divide by 10000 to convert to hectares).
  f -= log(square(80));///10000);
 
  // Pr(capture) calculated above.
  f -= log(pm + DBL_MIN);

SEPARABLE_FUNCTION void g_cluster(int i, const dvariable& D, const dvariable& g0, const dvariable& sigma, const dvariable& Xi, const dvariable& Yi)
  // Setting up variables
  int j,k;
  dvariable d,p,lambda,L1,L2,L3,indp,intL;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  // if (i==1){
  //   // Probabilities of caputure at each location for each trap.
  //   // Add a small amount to prevent zeros.
  //   p1=g0*mfexp(-square(dist)*pow(2*square(sigma),-1))+DBL_MIN; // Why doesn't -square(dist)/(2*square(sigma)) work?!
  //   p2=1-p1;
  //   // Probability of detection at any trap for each location.
  //   // Used for numerical integration of *detectable* animal density.
  //   for(j=1; j<=nmask; j++){
  //     p=1;
  //     for(k=1; k<=ntraps; k++){
  //       p*=p2(k)(j);
  //     }
  //     pm(j)=1-p;
  //   }
  //   lambda=A*D*sum(pm);
  //   //L2=-n*log(D*sum(pm));
  //   //L3=log_density_poisson(n,lambda);
  //   f -= -lambda;
  // }
  // Probability of capture histories for each animal.
  wi1=row(capt,i);
  // Prior on animal location.
  L1=log(D);
  for(j=1; j<=ntraps; j++){
    // Calculating distance between animal i and trap j.
    d=sqrt(square(Xi-trapsX(j))+square(Yi-trapsY(j)) + DBL_MIN);
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
  #include "minfil_mod.cpp"
  #include "df1b2gh.cpp"
