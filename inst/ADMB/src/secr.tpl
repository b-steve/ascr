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
  // Declaring indices.
  int i
  int j
  int u
  int k
  number dist
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
  number sum_probs
  number undet_prob
  number capt_prob
  matrix log_capt_probs(1,n_traps,1,n_mask)
  matrix log_evade_probs(1,n_traps,1,n_mask)
  matrix expected_ss(1,nr_expected_ss,1,nc_expected_ss)
  number ss_resid
  vector capt_hist(1,n_traps)

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
  // Start of likelihood calculation.
  f = 0.0;
  // Calculating expected signal strengths.
  if (fit_ss){
    if (linkfn_id == 1){
      expected_ss = detpars(1) - detpars(2)*dists;
    } else if (linkfn_id == 2){
      expected_ss = mfexp(detpars(1) - detpars(2)*dists);
    } else {
      cerr << "linkfn_id not recognised." << endl;
    }
  }
  // Calculating mask detection probabilities.
  sum_probs = 0;
  for (i = 1; i <= n_mask; i++){
    undet_prob = 1;
    for (j = 1; j <= n_traps; j++){
      dist = dists(j, i);
      if (fit_ss){
        ss_resid = cutoff - expected_ss(j, i);
      } else {
        ss_resid = 0;
      }
      capt_prob = detfn(dist, detpars, ss_resid);
      // Compare to calculating these outside loop.
      log_capt_probs(j, i) = log(capt_prob + DBL_MIN);
      log_evade_probs(j, i) = log(1 - capt_prob + DBL_MIN);
      undet_prob *= 1 - capt_prob;
    }
    sum_probs += 1 - undet_prob + DBL_MIN;
  }
  // Resetting i index.
  i = 0;
  // Contribution due to capture history.
  for (u = 1; u <= n_unique; u++){
    capt_hist = row(capt_bin_unique, u);
    if (local == 1){
      n_local = all_n_local(u);
      nc_localmats = n_local;
    }
    // Getting capture probabilities at local mask points.
    dvar_matrix * log_capt_probs_pointer;
    dvar_matrix * log_evade_probs_pointer;
    dvar_matrix local_log_capt_probs(1,nr_localmats,1,nc_localmats);
    dvar_matrix local_log_evade_probs(1,nr_localmats,1,nc_localmats);
    if (local){
      for (j = 1; j <= n_local; j++){
        local_log_capt_probs.colfill(j, column(log_capt_probs, all_which_local(u, j)));
        local_log_evade_probs.colfill(j, column(log_evade_probs, all_which_local(u, j)));
      } 
      log_capt_probs_pointer = &local_log_capt_probs;
      log_evade_probs_pointer = &local_log_evade_probs;
    } else {
      log_capt_probs_pointer = &log_capt_probs;
      log_evade_probs_pointer = &log_evade_probs;
    }
    dvar_vector bincapt_contrib(1,n_local);
    dvar_vector evade_contrib(1,n_local);
    evade_contrib = (1 - capt_hist)*(*log_evade_probs_pointer);
    // Calculating contribution due to uth unique capture history.
    if (fit_ss){      
      nr_local_expected_ss = nr_localmats;
      nc_local_expected_ss = nc_localmats;
    } else {
      bincapt_contrib = capt_hist*(*log_capt_probs_pointer) + evade_contrib;
      nr_local_expected_ss = 1;
      nc_local_expected_ss = 1;
    }
    // Filling local_expected_ss.
    dvar_matrix * expected_ss_pointer;
    dvar_matrix local_expected_ss(1,nr_local_expected_ss,1,nc_local_expected_ss);
    if (fit_ss){
      if (local){
        for (j = 1; j <= n_local; j++){
          local_expected_ss.colfill(j, column(expected_ss, all_which_local(u, j)));
        }
	expected_ss_pointer = &local_expected_ss;
      } else {
	expected_ss_pointer = &expected_ss;
      }
    }
    // Filling local_angs.
    // Move this sort of rubbish to the DATA_SECTION.
    if (fit_angs == 1){
      nr_local_angmat = nr_localmats;
      nc_local_angmat = nc_localmats;
    } else {
      nr_local_angmat = 1;
      nc_local_angmat = 1;
    }
    dmatrix * angs_pointer;
    dmatrix local_angs(1,nr_local_angmat,1,nc_local_angmat);
    if (fit_angs == 1){
      if (local){
        for (j = 1; j <= n_local; j++){
          local_angs.colfill(j, column(angs, all_which_local(u, j)));
        }
	angs_pointer = &local_angs;
      } else {
        angs_pointer = &angs;
      }
    }
    // Filling local_dists.
    if (fit_dists == 1){
      nr_local_dist = nr_localmats;
      nc_local_dist = nc_localmats;
    } else {
      nr_local_dist = 1;
      nc_local_dist = 1;
    }
    dmatrix * dists_pointer;
    dmatrix local_dists(1,nr_local_dist,1,nc_local_dist);
    if (fit_dists == 1){
      if (local){
        for (j = 1; j <= n_local; j++){
          local_dists.colfill(j, column(dists, all_which_local(u, j)));
        }
        dists_pointer = &local_dists;
      } else {
        dists_pointer = &dists;
      }
    }
    // Filling local_toa_ssq.
    if (fit_toas){
      nr_local_toa_ssq = n;
      nc_local_toa_ssq = nc_localmats;
    } else {
      nr_local_toa_ssq = 1;
      nc_local_toa_ssq = 1;
    }
    dmatrix * toa_ssq_pointer;
    dmatrix local_toa_ssq(1,nr_local_toa_ssq,1,nc_local_toa_ssq);
    if (fit_toas){
      if (local){
        for (j = 1; j <= n_local; j++){
          local_toa_ssq.colfill(j, column(toa_ssq, all_which_local(u, j)));
        }
        toa_ssq_pointer = &local_toa_ssq;
      } else {
        toa_ssq_pointer = &toa_ssq;
      }
    }
    for (k = 1; k <= capt_bin_freqs(u); k++){
      i++;
      // Contribution from capture locations.
      if (fit_ss){
        dvar_matrix local_log_ss_density(1,n_traps,1,n_local);
        for (j = 1; j <= n_traps; j++){
          local_log_ss_density.rowfill(j, log_dnorm(capt_ss(i, j), row((*expected_ss_pointer), j), detpars(3)));
        }
        bincapt_contrib = capt_hist*local_log_ss_density + evade_contrib;
      }
      if (any_suppars){
        dvar_vector supp_contrib(1,n_local);
        supp_contrib = 0;
        for (j = 1; j <= n_traps; j++){
          //  Try setting up a ragged array of capture locations for each individual instead.
          if (capt_bin_unique(u, j)){
            // Contribution from bearings.
            if (fit_angs){
              supp_contrib += suppars(kappa_ind)*cos(capt_ang(i, j) - row((*angs_pointer), j));
            }
            // Contribution from distances.
            if (fit_dists){
	      // Try saving alpha separately.
	      dvar_vector beta(1, n_local);
	      beta = suppars(alpha_ind)/row((*dists_pointer), j);
              supp_contrib += suppars(alpha_ind)*log(beta) + (suppars(alpha_ind) - 1)*log(capt_dist(i, j)) - beta*capt_dist(i, j);
	    }
          }
        }
        // Constant part of von Mises density.
        if (fit_angs){
          supp_contrib -= sum(capt_hist)*log(2*M_PI*bessi0(suppars(kappa_ind)));
        }
        // Constant part of the gamma density.
        if (fit_dists){
          supp_contrib -= sum(capt_hist)*gammln(suppars(alpha_ind));
        }
        // Contribution from times of arrival.
        if (fit_toas){
          // Try saving sigma_toa separately.
          supp_contrib += (1 - sum(capt_hist))*log(suppars(sigma_toa_ind)) - (row((*toa_ssq_pointer), i)/(2*square(suppars(sigma_toa_ind))));
        }
        f -= log(sum(mfexp(bincapt_contrib + supp_contrib)) + DBL_MIN);
      } else {
        f -= log(sum(mfexp(bincapt_contrib)) + DBL_MIN);
      }
    }
  }
  // Calculating ESA.
  esa = A*sum_probs;
  // Contribution from n.
  f -= log_dpois(n, D*esa);
  // Extra bit that falls out of ll.
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
  #include <detfuns.cpp>
  #include <invlinkfuns.cpp>

