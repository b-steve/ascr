DATA_SECTION
  init_int n
  init_int ntraps
  init_int nmask
  init_number A
  init_number trace
  init_vector D_bounds(1,2)
  init_number D_phase
  init_number n_detpars
  init_vector detpars_lb(1,n_detpars)
  init_vector detpars_ub(1,n_detpars)
  init_vector detpars_phase(1,n_detpars)
  init_vector suppars_lb(1,3)
  init_vector suppars_ub(1,3)
  init_vector suppars_phase(1,3)
  init_number detfn_id
  init_number DBL_MIN
  init_matrix capt(1,n,1,ntraps)
  init_matrix dist(1,ntraps,1,nmask)

PARAMETER_SECTION
  objective_function_value f
  init_bounded_number D(D_bounds[1],D_bounds[2],D_phase)
  init_bounded_number_vector detpars(detpars_lb,detpars_ub,detpars_phase)
  init_bounded_number_vector suppars(suppars_lb,suppars_ub,suppars_phase)

PROCEDURE_SECTION
  const double pi = 3.141592653589793238463;
  f = 0.0;  

GLOBALS_SECTION
  #include <densfuns.cpp>
  #include <detfuns.cpp>