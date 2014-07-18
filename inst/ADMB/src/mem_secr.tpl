DATA_SECTION
  init_int n_unique
  init_int local
  // Ragged matrix of local point indices.
  init_ivector all_n_local(1,n_unique)
  init_matrix all_which_local(1,n_unique,1,all_n_local)
  // Density parameter details.
  init_number D_lb
  init_number D_ub
  init_number D_phase
  init_number D_sf
  // Detection function parameter details.
  init_int n_detpars
  init_vector detpars_lb(1,n_detpars)
  init_vector detpars_ub(1,n_detpars)
  init_ivector detpars_phase(1,n_detpars)
  init_vector detpars_sf(1,n_detpars)
  init_ivector detpars_linkfns(1,n_detpars)
  // Supplementary info parameter details.
  init_int n_suppars
  init_vector suppars_lb(1,n_suppars)
  init_vector suppars_ub(1,n_suppars)
  init_ivector suppars_phase(1,n_suppars)
  init_vector suppars_sf(1,n_suppars)
  init_ivector suppars_linkfns(1,n_suppars)
  // Calculating total number of estimated parameters.
  int n_ests
  !! n_ests = 0;
  !! if (D_phase > -1){
  !!   n_ests++;
  !! }
  !! for (i = 1; i <= n_detpars; i++){
  !!   if (detpars_phase(i) > -1){
  !!     n_ests++;
  !!   }
  !! }
  !! for (i = 1; i <= n_suppars; i++){
  !!   if (suppars_phase(i) > -1){
  !!     n_ests++;
  !!   }
  !! }
  // Initialising various scalars.
  init_number detfn_id
  init_number buffer
  init_number trace
  init_number DBL_MIN
  init_int n
  init_int n_traps
  init_int n_mask
  init_number A
  // Initialising capture histories and frequencies.
  init_matrix capt_bin_unique(1,n_unique,1,n_traps)
  init_ivector capt_bin_freqs(1,n_unique)
  // Initialising observed bearings.
  init_int fit_angs
  int nr_ang
  int nc_ang
  int nr_angmat
  int nc_angmat
  int nr_local_angmat
  int nc_local_angmat
  !! if (fit_angs == 1){
  !!   nr_ang = n;
  !!   nc_ang = n_traps;
  !!   nr_angmat = n_traps;
  !!   nc_angmat = n_mask;
  !! } else {
  !!   nr_ang = 1;
  !!   nc_ang = 1;
  !!   nr_angmat = 1;
  !!   nc_angmat = 1;
  !! }
  init_matrix capt_ang(1,nr_ang,1,nc_ang)
  // Initialising observed distances.
  init_int fit_dists
  int nr_dist
  int nc_dist
  int nr_local_dist
  int nc_local_dist
  !! if (fit_dists == 1){
  !!   nr_dist = n;
  !!   nc_dist = n_traps;
  !! } else {
  !!   nr_dist = 1;
  !!   nc_dist = 1;
  !! }
  init_matrix capt_dist(1,nr_dist,1,nc_dist)
  // Initialising observed signal strengths.
  init_int fit_ss
  init_number cutoff
  init_int linkfn_id
  int nr_ss
  int nc_ss
  int nr_expected_ss
  int nc_expected_ss
  int nr_local_expected_ss
  int nc_local_expected_ss
  !! if (fit_ss == 1){
  !!   nr_ss = n;
  !!   nc_ss = n_traps;
  !!   nr_expected_ss = n_traps;
  !!   nc_expected_ss = n_mask;
  !! } else {
  !!   nr_ss = 1;
  !!   nc_ss = 1;
  !!   nr_expected_ss = 1;
  !!   nc_expected_ss = 1;
  !! }
  init_matrix capt_ss(1,nr_ss,1,nc_ss)
  // Initialising observed times of arrival.
  init_int fit_toas
  int nr_toa
  int nc_toa
  int nr_toa_ssq
  int nc_toa_ssq
  int nr_local_toa_ssq
  int nc_local_toa_ssq
  !! if (fit_toas == 1){
  !!   nr_toa = n;
  !!   nc_toa = n_traps;
  !!   nr_toa_ssq = n;
  !!   nc_toa_ssq = n_mask;
  !! } else {
  !!   nr_toa = 1;
  !!   nc_toa = 1;
  !!   nr_toa_ssq = 1;
  !!   nc_toa_ssq = 1;
  !! }
  init_matrix capt_toa(1,nr_toa,1,nc_toa)
  // Initialising observed exact distances (not yet implemented?).
  init_int fit_mrds
  int nr_mrds
  int nc_mrds
  !! if (fit_mrds == 1){
  !!   nr_mrds = n;
  !!   nc_mrds = n_traps;
  !! } else {
  !!   nr_mrds = 1;
  !!   nc_mrds = 1;
  !! }
  // Expected supplementary information from each mask point.
  init_matrix mrds_dist(1,nr_mrds,1,nc_mrds)
  init_matrix dists(1,n_traps,1,n_mask)
  init_matrix angs(1,nr_angmat,1,nc_angmat)
  init_matrix toa_ssq(1,nr_toa_ssq,1,nc_toa_ssq)
  // Indicator for supplementary information (other than SS).
  int any_suppars
  !! if (fit_angs + fit_dists + fit_toas > 0){
  !!   any_suppars = 1;
  !! } else {
  !!   any_suppars = 0;
  !! }
  // Sorting out positions in suppars.
  int kappa_ind
  int alpha_ind
  int sigma_toa_ind
  int curr_ind
  !! curr_ind = 1;
  !! if (fit_angs){
  !!   kappa_ind = curr_ind;
  !!   curr_ind++;
  !! }
  !! if (fit_dists){
  !!   alpha_ind = curr_ind;
  !!   curr_ind++;
  !! }
  !! if (fit_toas){
  !!   sigma_toa_ind = curr_ind;
  !! }
  // Setting n_local permanently if local integration is disabled.
  int nr_localmats
  int nc_localmats
  int n_local
  !! if (local == 0){
  !!   n_local = n_mask;
  !!   nr_localmats = 1;
  !!   nc_localmats = 1;
  !! } else {
  !!   nr_localmats = n_traps;
  !! }
  // Non-differentiable objects for PROCEDURE_SECTION.
  int i
  int j
  int k
  int m
  int t
  number dist
  number n_dets
  int n_u_contribs
  !! if (fit_ss){
  !!   n_u_contribs = n;
  !! } else {
  !!   n_u_contribs = n_unique;
  !! }

PARAMETER_SECTION
  objective_function_value f
  init_bounded_number D_link(D_lb,D_ub,D_phase)
  init_bounded_number_vector detpars_link(1,n_detpars,detpars_lb,detpars_ub,detpars_phase)
  init_bounded_number_vector suppars_link(1,n_suppars,suppars_lb,suppars_ub,suppars_phase)
  !! D_link.set_scalefactor(D_sf);
  !! detpars_link.set_scalefactor(detpars_sf);
  !! suppars_link.set_scalefactor(suppars_sf);
  sdreport_vector par_ests(1,n_ests)
  sdreport_number esa
  number D
  vector detpars(1,n_detpars)
  vector suppars(1,n_suppars)
  vector log_u_contribs(1,n_u_contribs)
  vector i_contribs(1,n)
  number log_s_contribs
  number log_b_contribs
  number point_capt
  number point_evade
  number sum_probs
  number capt_prob

PROCEDURE_SECTION
  // Grabbing detection function.
  detfn_pointer detfn = get_detfn(detfn_id);
  // Converting linked parameters to real parameters and setting up par_ests vector.
  D = mfexp(D_link);
  j = 1;
  if (D_phase > -1){
    par_ests(j) = D;
    j++;
  }
  // Converting parameters from their link scales.
  invlinkfn_pointer invlinkfn;
  for (i = 1; i <= n_detpars; i++){
    invlinkfn = get_invlinkfn(detpars_linkfns(i));
    detpars(i) = invlinkfn(detpars_link(i));
    if (detpars_phase(i) > -1){
      par_ests(j) = detpars(i);
      j++;
    }
  }
  for (i = 1; i <= n_suppars; i++){
    invlinkfn = get_invlinkfn(suppars_linkfns(i));
    suppars(i) = invlinkfn(suppars_link(i));
    if (suppars_phase(i) > -1){
      par_ests(j) = suppars(i);
      j++;
    }
  }
  // Looping over mask points. Simultaneously calculating
  // contributions from each detection.
  i_contribs = 0;
  sum_probs = 0;
  for (m = 1; m <= n_mask; m++){
    log_u_contribs = 0;
    point_evade = 1;
    // Calculating contribution due to capture location.
    for (t = 1; t <= n_traps; t++){
      dist = dists(t, m);
      capt_prob = detfn(dist, detpars, cutoff);
      if (fit_ss){
        dvariable expected_ss;
        expected_ss = detpars(1) - detpars(2)*dist;
        k = 0;
        for (i = 1; i <= n_unique; i++){
          for (j = 1; j <= capt_bin_freqs(i); j++){
            k++;
            if (capt_bin_unique(i, t)){
              log_u_contribs(k) += log_dnorm(capt_ss(k, t), expected_ss, detpars(3));
            } else {
              log_u_contribs(k) += log(1 - capt_prob + DBL_MIN);
            }
          }
        }
      } else {
        for (i = 1; i <= n_unique; i++){
          if (capt_bin_unique(i, t) == 1){
            log_u_contribs(i) += log(capt_prob + DBL_MIN);
          } else {
            log_u_contribs(i) += log(1 - capt_prob + DBL_MIN);
          }
        }
      }
      point_evade *= 1 - capt_prob;
    }
    point_capt = 1 - point_evade;
    sum_probs += point_capt;
    k = 0;
    for (i = 1; i <= n_unique; i++){
      n_dets = sum(row(capt_bin_unique, i));
      for (j = 1; j <= capt_bin_freqs(i); j++){
        k++;
        log_s_contribs = 0;
        if (fit_toas){
          if (n_dets > 1){
            log_s_contribs += (1 - n_dets)*log(suppars(sigma_toa_ind)) - toa_ssq(k, m)/(2*square(suppars(sigma_toa_ind)));
          }
        }
        if (fit_ss){
          log_b_contribs = log_u_contribs(k);
        } else {
          log_b_contribs = log_u_contribs(i);
        }
        i_contribs(k) += mfexp(log_b_contribs + log_s_contribs);
      }
    }
  }
  esa = A*sum_probs;
  // Contribution from capture histories.
  f = -sum(log(i_contribs + DBL_MIN));
  // Contribution from n.
  f -= log_dpois(n, D*esa);
  // Extra bit that falls out of log-likelihood.
  f -= -n*log(sum_probs);
  // Printing trace.
  if (trace){
    cout << "D: " << D << ", ";
    for (i = 1; i <= n_detpars; i++){
      cout << "DF Par " << i << ": " << detpars(i) << ", ";
    }
    if (any_suppars){
      for (i = 1; i <= n_suppars; i++){
        cout << "Supp Par " << i << ": " << suppars(i) << ", ";
      }
    }
    cout << "LL: " << -f << endl;
  }

GLOBALS_SECTION
  #include <densfuns.cpp>
  #include <mem_detfuns.cpp>
  #include <invlinkfuns.cpp>
