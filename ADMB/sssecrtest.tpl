TOP_OF_MAIN_SECTION
  arrmblsize=10000000;

DATA_SECTION

  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_matrix sscapt(1,345,1,6)
  init_matrix dist(1,6,1,3449)
  init_matrix capt(1,345,1,6)

PARAMETER_SECTION

  objective_function_value f
  init_bounded_number D(0,1e+07)
  init_bounded_number g0(0,1)
  init_bounded_number sigma(0,1e+07)
  init_number ssb0
  init_number ssb1
  init_bounded_number sigmass(0,6e+05)

PROCEDURE_SECTION
  // Setting up variables
  const double pi=3.14159265359;
  const double c=150;
  int i,j,k;
  //dvariable p,lambda,L1,L2,L3;
  dvariable muss,p,lambda,L1,L2,L3;
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector ci1(1,ntraps);
  dvar_vector ssll(1,nmask);
  for(i=1; i<=nmask; i++){
    p=1;
    for(j=1; j<=ntraps; j++){
      muss=ssb0+ssb1*dist(j,i);
      p*=cumd_norm((log(c)-muss)/sigmass);
      //p*=1-(g0*mfexp(-square(dist(j,i))/(2*square(sigma)))+DBL_MIN);
    }
    pm(i)=1-p;
  }
  L1=0;
  // Probability of capture histories for each animal.
  for(i=1; i<=n; i++){
  ssll=0;
    wi1=row(sscapt,i);
    ci1=row(capt,i);
    //wi1=row(capt,i);
    for(j=1; j<=nmask; j++){
      p=1;
      for(k=1; k<=ntraps; k++){
        muss=ssb0+ssb1*dist(k,j);
        if(ci1(k)==0){
          //p*=1-(g0*mfexp(-square(dist(k,j))/(2*square(sigma)))+DBL_MIN);
          p*=cumd_norm((log(c)-muss)/sigmass);
        }else{
          //p*=g0*mfexp(-square(dist(k,j))/(2*square(sigma)))+DBL_MIN;
          //p*=1-cumd_norm((log(c)-muss)/sigmass);
          p*=mfexp(square(log(wi1(k))-muss)/(2*square(sigmass)))/(sigmass*sqrt(2*pi));
        }
        //p*=1-cumd_norm((log(c)-muss)/sigmass);
        //p*=g0*mfexp(-square(dist(k,j))/(2*square(sigma)))+DBL_MIN;
      }
      ssll(j)=D*p;
    }
    L1+=log(sum(ssll)+DBL_MIN);
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm)+DBL_MIN;
  L2=-n*log(D*sum(pm)+DBL_MIN);
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  cout << L1 << " " << L2 << " " << L3 << endl;

GLOBALS_SECTION
  #include <float.h>

REPORT_SECTION
