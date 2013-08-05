  #include <float.h>
#include <admodel.h>
#include <contrib.h>

#include <df1b2fun.h>

#include <adrndeff.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <sepsssecrw.htp>

  df1b2_parameters * df1b2_parameters::df1b2_parameters_ptr=0;
  model_parameters * model_parameters::model_parameters_ptr=0;
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
  model_parameters_ptr=this;
  initializationfunction();
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
  D.allocate(0,1e+07,"D");
  ssb0.allocate("ssb0");
  ssb1.allocate(-10,0,"ssb1");
  sigmass.allocate(0,1e+07,"sigmass");
  X.allocate(1,n,"X");
  Y.allocate(1,n,"Y");
}
void model_parameters::userfunction(void)
{
  f =0.0;
  int i;
  f=0;
  for(i=1;i<=n;i++){
    g_cluster(i, D, ssb0, ssb1, sigmass, X(i), Y(i));
  }
}

void SEPFUN1  model_parameters::g_cluster(int i, const dvariable& D, const dvariable& ssb0, const dvariable& ssb1, const dvariable& sigmass, const dvariable& x_i, const dvariable& y_i)
{
  begin_df1b2_funnel();
  f -= -0.5*square(x_i);
  f -= -0.5*square(y_i);
  // Varialbes with uniform distribution on [-40,40]
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
    f -= -lambda;
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
  end_df1b2_funnel();
}
  long int arrmblsize=0;

int main(int argc,char * argv[])
{
  ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize=150000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_MAX_NVAR_OFFSET(4600404);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    df1b2variable::noallocate=1;
    df1b2_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;

    function_minimizer::random_effects_flag=1;
    df1b2variable::noallocate=0;
    mp.preliminary_calculations();
    initial_df1b2params::separable_flag=1;
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

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
  }

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(void){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

void df1b2_parameters::user_function(void)
{
  f =0.0;
  int i;
  f=0;
  for(i=1;i<=n;i++){
    g_cluster(i, D, ssb0, ssb1, sigmass, X(i), Y(i));
  }
}

void   df1b2_pre_parameters::g_cluster(int i, const funnel_init_df1b2variable& D, const funnel_init_df1b2variable& ssb0, const funnel_init_df1b2variable& ssb1, const funnel_init_df1b2variable& sigmass, const funnel_init_df1b2variable& x_i, const funnel_init_df1b2variable& y_i)
{
  begin_df1b2_funnel();
  f -= -0.5*square(x_i);
  f -= -0.5*square(y_i);
  // Varialbes with uniform distribution on [-40,40]
  df1b2variable Xi = 80.0*cumd_norm(x_i) - 40.0;
  df1b2variable Yi = 80.0*cumd_norm(y_i) - 40.0;
  
  // Setting up variables
  int j,k;
  const double pi=3.14159265359;
  df1b2variable d,p,lambda,L1,L2,L3,exprss;
  df1b2matrix p1(1,ntraps,1,nmask);
  df1b2matrix p2(1,ntraps,1,nmask);
  df1b2matrix muss(1,ntraps,1,nmask);
  df1b2vector pm(1,nmask);
  df1b2vector wi1(1,ntraps);
  df1b2vector ci1(1,ntraps);
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
    f -= -lambda;
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
  end_df1b2_funnel();
}
   
void df1b2_pre_parameters::setup_quadprior_calcs(void) 
{ 
  df1b2_gradlist::set_no_derivatives(); 
  quadratic_prior::in_qp_calculations=1; 
}  
  
void df1b2_pre_parameters::begin_df1b2_funnel(void) 
{ 
  (*re_objective_function_value::pobjfun)=0; 
  other_separable_stuff_begin(); 
  f1b2gradlist->reset();  
  if (!quadratic_prior::in_qp_calculations) 
  { 
    df1b2_gradlist::set_yes_derivatives();  
  } 
  funnel_init_var::allocate_all();  
}  
 
void df1b2_pre_parameters::end_df1b2_funnel(void) 
{  
  lapprox->do_separable_stuff(); 
  other_separable_stuff_end(); 
} 
  
void model_parameters::begin_df1b2_funnel(void) 
{ 
  if (lapprox)  
  {  
    {  
      begin_funnel_stuff();  
    }  
  }  
}  
 
void model_parameters::end_df1b2_funnel(void) 
{  
  if (lapprox)  
  {  
    end_df1b2_funnel_stuff();  
  }  
} 

void df1b2_parameters::allocate(void) 
{
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  f.allocate("f");  /* ADOBJECTIVEFUNCTION */
  D.allocate(0,1e+07,"D");
  ssb0.allocate("ssb0");
  ssb1.allocate(-10,0,"ssb1");
  sigmass.allocate(0,1e+07,"sigmass");
  X.allocate(1,n,"X");
  Y.allocate(1,n,"Y");
}
