  #include <float.h>
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <secr_gen.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  ntraps.allocate("ntraps");
  nmask.allocate("nmask");
  A.allocate("A");
  capt.allocate(1,345,1,6,"capt");
  dist.allocate(1,6,1,3449,"dist");
  trace.allocate("trace");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  D.allocate(0,1e+08,"D");
  par0.allocate(-50,50,"par0");
  par1.allocate(-10,0,"par1");
}

void model_parameters::userfunction(void)
{
  f =0.0;
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
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
