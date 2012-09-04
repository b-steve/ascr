TOP_OF_MAIN_SECTION
  arrmblsize=150000000;

DATA_SECTION

  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_matrix capt(1,345,1,6)
  init_matrix dist(1,6,1,3449)
  init_matrix traps(1,6,1,2)
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
  f=0;
  cout << D << " " << g0 << " " << sigma << " " << f << endl;
  //cout << X << endl;
  for(int i=1; i<=n; i++){
    g_cluster(i, D, g0, sigma, X(i), Y(i));
  }

SEPARABLE_FUNCTION void g_cluster(int i, const dvariable& D, const dvariable& g0, const dvariable& sigma, const dvariable& Xi, const dvariable& Yi)
  // Setting up variables
  int j,k;
  dvariable p,lambda,L1,L2,L3,d,pc;
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
  p1=g0*mfexp(-square(dist)*pow(2*square(sigma),-1))+DBL_MIN;
  p2=1-p1;
  logp1=log(p1);
  logp2=log(p2);
  // Probability of detection at any trap for each location.
  for(j=1; j<=nmask; j++){
    p=1;
    for(k=1; k<=ntraps; k++){
      p*=p2(k)(j);
    }
    pm(j)=1-p;
  }
  // Probability of capture histories for each animal.
  wi1=capt(i)(1,ntraps);
  p=1;
  if (i == 1){
  //cout << trapsX << endl;
  //cout << trapsY << endl;
  }
  for (j=1; j<=ntraps; j++){
    d=pow(square(Xi-trapsX(j))+square(Yi-trapsY(j)),0.5);
    pc=g0*mfexp(-square(d)/(2*square(sigma)))+DBL_MIN;
    p*=wi1(j)*pc+(1-wi1(j))*(1-pc);
  }
  L1=log(D*p);
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2/n+L3/n);
  if (i == 1){
    //cout << p << endl;
    //cout << L1 << " " << L2 << " " << L3 << endl;
  }

GLOBALS_SECTION
  #include <float.h>

