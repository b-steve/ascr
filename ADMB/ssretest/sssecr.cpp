  #include <float.h>
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <sssecr.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  n.allocate("n");
  ntraps.allocate("ntraps");
  nmask.allocate("nmask");
  A.allocate("A");
  capt.allocate(1,40,1,6,"capt");
  sscapt.allocate(1,40,1,6,"sscapt");
  dist.allocate(1,6,1,3449,"dist");
  c.allocate("c");
  traps.allocate(1,6,1,2,"traps");
  trace.allocate("trace");
  trapsX.allocate(1,ntraps);
  trapsY.allocate(1,ntraps);
trapsX=column(traps,1);
trapsY=column(traps,2);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  D.allocate(0,1e+07,"D");
  ssb0.allocate("ssb0");
  ssb1.allocate(-10,0,"ssb1");
  sigmass.allocate(0,1e+07,"sigmass");
}

void model_parameters::userfunction(void)
{
  f =0.0;
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
