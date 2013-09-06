DATA_SECTION
  init_number D_lb
  init_number D_ub
  init_number D_phase
  init_number D_sf
  init_int n_detpars
  init_vector detpars_lb(1,n_detpars)
  init_vector detpars_ub(1,n_detpars)
  init_ivector detpars_phase(1,n_detpars)
  init_vector detpars_sf(1,n_detpars)
  init_int n_suppars
  init_vector suppars_lb(1,3)
  init_vector suppars_ub(1,3)
  init_ivector suppars_phase(1,3)
  init_vector suppars_sf(1,n_suppars)
  init_number detfn_id
  init_number trace
  init_number DBL_MIN
  init_int n
  init_int n_traps
  init_int n_mask
  init_number A
  init_number n_freqs
  init_vector call_freqs(1,n_freqs)
  init_matrix capt_bin(1,n,1,n_traps)
  init_ivector ang_dim(1,2)
  int nr_ang
  int nc_ang
  !! nr_ang = ang_dim[1];
  !! nc_ang = ang_dim[2];
  init_matrix capt_ang(1,nr_ang,1,nc_ang);
  init_ivector dist_dim(1,2)
  int nr_dist
  int nc_dist
  !! nr_dist = dist_dim[1];
  !! nc_dist = dist_dim[2];
  init_matrix capt_dist(1,nr_dist,1,nc_dist);
  init_ivector ss_dim(1,2)
  int nr_ss
  int nc_ss
  !! nr_ss = ss_dim[1];
  !! nc_ss = ss_dim[2];
  init_matrix capt_ss(1,nr_ss,1,nc_ss);
  init_ivector toa_dim(1,2)
  int nr_toa
  int nc_toa
  !! nr_toa = toa_dim[1];
  !! nc_toa = toa_dim[2];
  init_matrix capt_toa(1,nr_toa,1,nc_toa)
  init_ivector mrds_dim(1,2)
  int nr_mrds
  int nc_mrds
  !! nr_mrds = mrds_dim[1];
  !! nc_mrds = mrds_dim[2];
  init_matrix mrds_dist(1,nr_mrds,1,nc_mrds);
  init_matrix dists(1,n_traps,1,n_mask)

PARAMETER_SECTION
  objective_function_value f
  init_bounded_number D(D_lb,D_ub,D_phase)
  init_bounded_number_vector detpars(1,n_detpars,detpars_lb,detpars_ub,detpars_phase)
  init_bounded_number_vector suppars(1,n_suppars,suppars_lb,suppars_ub,suppars_phase)
  !! D.set_scalefactor(D_sf);
  !! detpars.set_scalefactor(detpars_sf);
  !! suppars.set_scalefactor(suppars_sf);

PROCEDURE_SECTION
  detfn_pointer detfn = get_detfn(detfn_id);
  const double pi = 3.141592653589793238463;
  f = 0.0;
  // Calculating mask detection probabilities.
  double dist;
  dvariable sum_probs, undet_prob, capt_prob;
  dvar_matrix log_capt_probs(1, n_traps, 1, n_mask);
  dvar_matrix log_evade_probs(1, n_traps, 1, n_mask);
  sum_probs = 0;
  for (int i = 1; i <= n_mask; i++){
    undet_prob = 1;
    for (int j = 1; j <= n_traps; j++){
      dist = dists(j, i);
      capt_prob = detfn(dist, detpars);
      // Compare to calculating these outside loop.
      log_capt_probs[j, i] = log(capt_prob + DBL_MIN);
      log_evade_probs[j, i] = log(1 - capt_prob + DBL_MIN);
      undet_prob *= 1 - capt_prob;
    }
    sum_probs += 1 - undet_prob;
  }
  dvar_vector capt_hist(1, n_traps);
  dvar_vector secr_contrib(1, n_mask);
  for (int i = 1; i <= n_mask; i++){
    f -= log(sum(mfexp(capt_hist*log_capt_probs + (1 - capt_hist)*log_evade_probs)) + DBL_MIN);
  }
  // Contribution from n.
  f -= log_dpois(n, A*D*sum_probs);
  f -= -n*log(D*sum_probs);
  // Printing trace.
  if (trace){
    cout << "D: " << D << ", ";
    for (int i; i <= n_detpars; i++){
      cout << "DF Par " << i << ": " << detpars[i] << ", ";
    }
    for (int i; i <= n_suppars; i++){
      cout << "Supp Par " << i << ": " << detpars[i] << ", ";
    }
    cout << "LL: " << -f << endl;
  }

GLOBALS_SECTION
  #include <densfuns.cpp>
  #include <detfuns.cpp>

