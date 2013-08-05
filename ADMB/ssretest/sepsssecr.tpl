// Separable function for a(theta).
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
  init_number ssb0(-1)
  init_bounded_number ssb1(-10,0,-1)
  init_bounded_number sigmass(1,1e+07,-1)
  random_effects_vector X(1,n)
  random_effects_vector Y(1,n)
  random_effects_vector m(1,2)

PROCEDURE_SECTION
  int i;
  f=0;
  ll_exp(D, ssb0, ssb1, sigmass, m(1), m(2));
  for(i=1;i<=n;i++){
    ll_ind(i, D, ssb0, ssb1, sigmass, X(i), Y(i));
  }

  // Exponential part of the likelihood.
SEPARABLE_FUNCTION void ll_exp(const dvariable& D, const dvariable& ssb0, const dvariable& ssb1, const dvariable& sigmass, const dvariable& m1, const dvariable& m2)
  const double pi = 3.14159265359;
  dvariable d, muss;
  f -= -0.5*square(m1) - 0.5*log(2*pi);
  f -= -0.5*square(m2) - 0.5*log(2*pi);
  // Variables with uniform distribution on [-40, 40]
  dvariable x = 80.0*cumd_norm(m1) - 40.0;
  dvariable y = 80.0*cumd_norm(m2) - 40.0;
  dvariable p = 1.0;
  for (int j = 1; j <= ntraps; j++){
    d = sqrt(square(x - trapsX(j)) + square(y - trapsY(j)));
    muss = ssb0 + ssb1*d;
    // Probability of evasion at trap j.
    p *= cumd_norm((c - muss)/sigmass);
  }
  // Probability of capture.
  dvariable pm = 1 - p;

  // Contribution to the likelihood is: exp(-1 * \int D * Area * Pr(capture) * f(X) dX)
  // Here we need to give: log(D * Area * Pr(capture) * f(X))
  f -= log(D);

  // Area of support of X = 80^2/10000 (divide by 10000 to convert to hectares).
  f -= log(square(80)/10000);
 
  // Pr(capture) calculated above.
  f -= log(pm + DBL_MIN);

SEPARABLE_FUNCTION void ll_ind(int i, const dvariable& D, const dvariable& ssb0, const dvariable& ssb1, const dvariable& sigmass, const dvariable& x_i, const dvariable& y_i)

  f -= -0.5*square(x_i);
  f -= -0.5*square(y_i);

  // Variables with uniform distribution on [-40,40]
  dvariable Xi = 80.0*cumd_norm(x_i) - 40.0;
  dvariable Yi = 80.0*cumd_norm(y_i) - 40.0;
  
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

  //Probability of capture histories for each animal.
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
      L1 += log(cumd_norm((c-exprss)/sigmass)+DBL_MIN);
    }
  }
  f-=L1;
  //cout << "L1: " << L1 << endl;

GLOBALS_SECTION
  #include <float.h>
  #include "minfil_mod.cpp"
  #include "df1b2ghmult.cpp"
  #include "df1b2f15.cpp"

TOP_OF_MAIN_SECTION
  arrmblsize=150000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_MAX_NVAR_OFFSET(4600404);
