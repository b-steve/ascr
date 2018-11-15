DATA_SECTION
  // Session-related stuff.
  init_int n_sessions
  init_vector survey_length(1,n_sessions)
  init_ivector n_unique_per_sess(1,n_sessions)
  init_int local
  // Ragged matrix of local point indices.
  init_imatrix n_local_per_unique(1,n_sessions,1,n_unique_per_sess)
  init_3darray which_local_per_unique(1,n_sessions,1,n_unique_per_sess,1,n_local_per_unique)
  // Density parameter details.
  init_int n_D_betapars
  init_vector D_betapars_lb(1,n_D_betapars)
  init_vector D_betapars_ub(1,n_D_betapars)
  init_ivector D_betapars_phase(1,n_D_betapars)
  init_vector D_betapars_sf(1,n_D_betapars)
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
  int b
  int a
  int c
  int d
  int s
  int i_start
  int i_end
  number dist
  // Calculating total number of estimated parameters.
  int n_ests
  !! n_ests = 0;
  !! for (i = 1; i <= n_D_betapars; i++){
  !!   if (D_betapars_phase(i) > -1){
  !!     n_ests++;
  !!   }
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
  init_number trace
  init_number dbl_min
  init_ivector n_per_sess(1,n_sessions)
  init_ivector n_traps_per_sess(1,n_sessions)
  init_ivector n_mask_per_sess(1,n_sessions)
  int n_mask_total
  !! n_mask_total = sum(n_mask_per_sess);
  init_vector A_per_sess(1,n_sessions)
  // Initialising capture histories and frequencies.
  init_3darray capt_bin_unique(1,n_sessions,1,n_unique_per_sess,1,n_traps_per_sess)
  init_imatrix capt_bin_freqs(1,n_sessions,1,n_unique_per_sess)
  // Initialising observed bearings.
  init_int fit_angs
  // Logical indicator for directional calling.
  init_int fit_dir
  init_number n_dir_quadpoints
  // Workaround for ADMB bug with 4D arrays.
  int dummy_n_dir_quadpoints
  !! dummy_n_dir_quadpoints = n_dir_quadpoints;
  init_int fit_het_source
  init_int het_source_gh
  init_number n_het_source_quadpoints
  int n_het_source_nodes
  !! if (fit_het_source & het_source_gh){
  !!   n_het_source_nodes = n_het_source_quadpoints;
  !! } else {
  !!   n_het_source_nodes = 1;
  !! }
  init_vector het_source_nodes(1,n_het_source_nodes)
  init_vector het_source_weights(1,n_het_source_nodes)
  ivector length_fs(1,n_sessions)
  !! if (fit_dir){
  !!   length_fs = n_per_sess;
  !! } else {
  !!   length_fs = 1;
  !! }
  ivector nr_ang(1,n_sessions)
  ivector nc_ang(1,n_sessions)
  ivector nr_angmat(1,n_sessions)
  ivector nc_angmat(1,n_sessions)
  int nr_local_angmat
  int nc_local_angmat
  !! if (fit_angs){
  !!   nr_ang = n_per_sess;
  !!   nc_ang = n_traps_per_sess;
  !!   nr_angmat = n_traps_per_sess;
  !!   nc_angmat = n_mask_per_sess;
  !! } else if (fit_dir) {
  !!   nr_ang = 1;
  !!   nc_ang = 1;
  !!   nr_angmat = n_traps_per_sess;
  !!   nc_angmat = n_mask_per_sess;
  !! } else {
  !!   nr_ang = 1;
  !!   nc_ang = 1;
  !!   nr_angmat = 1;
  !!   nc_angmat = 1;
  !! }
  init_3darray capt_ang(1,n_sessions,1,nr_ang,1,nc_ang)
  // Initialising observed distances.
  init_int fit_dists
  ivector nr_dist(1,n_sessions)
  ivector nc_dist(1,n_sessions)
  int nr_local_dist
  int nc_local_dist
  !! if (fit_dists == 1){
  !!   nr_dist = n_per_sess;
  !!   nc_dist = n_traps_per_sess;
  !! } else {
  !!   nr_dist = 1;
  !!   nc_dist = 1;
  !! }
  init_3darray capt_dist(1,n_sessions,1,nr_dist,1,nc_dist)
  // Initialising observed signal strengths.
  init_int fit_ss
  init_number cutoff
  init_int first_calls
  ivector n_mask_det_probs(1,n_sessions)
  !! if (first_calls == 1){
  !!   n_mask_det_probs = n_mask_per_sess;
  !! } else {
  !!   n_mask_det_probs = 1;
  !! }
  init_number lower_cutoff
  init_int linkfn_id
  ivector nr_ss(1,n_sessions)
  ivector nc_ss(1,n_sessions)
  ivector nr_expected_ss(1,n_sessions)
  ivector nc_expected_ss(1,n_sessions)
  ivector nr_local_expected_ss(1,n_sessions)
  ivector nc_local_expected_ss(1,n_sessions)
  !! if (fit_ss == 1){
  !!   nr_ss = n_per_sess;
  !!   nc_ss = n_traps_per_sess;
  !!   nr_expected_ss = n_traps_per_sess;
  !!   nc_expected_ss = n_mask_per_sess;
  !! } else {
  !!   nr_ss = 1;
  !!   nc_ss = 1;
  !!   nr_expected_ss = 1;
  !!   nc_expected_ss = 1;
  !! }
  init_3darray capt_ss(1,n_sessions,1,nr_ss,1,nc_ss)
  // Initialising observed times of arrival.
  init_int fit_toas
  ivector nr_toa(1,n_sessions)
  ivector nc_toa(1,n_sessions)
  ivector nr_toa_ssq(1,n_sessions)
  ivector nc_toa_ssq(1,n_sessions)
  ivector nr_local_toa_ssq(1,n_sessions)
  ivector nc_local_toa_ssq(1,n_sessions)
  !! if (fit_toas == 1){
  !!   nr_toa = n_per_sess;
  !!   nc_toa = n_traps_per_sess;
  !!   nr_toa_ssq = n_per_sess;
  !!   nc_toa_ssq = n_mask_per_sess;
  !! } else {
  !!   nr_toa = 1;
  !!   nc_toa = 1;
  !!   nr_toa_ssq = 1;
  !!   nc_toa_ssq = 1;
  !! }
  init_3darray capt_toa(1,n_sessions,1,nr_toa,1,nc_toa)
  // Initialising observed exact distances (not yet implemented?).
  init_int fit_mrds
  ivector nr_mrds(1,n_sessions)
  ivector nc_mrds(1,n_sessions)
  !! if (fit_mrds == 1){
  !!   nr_mrds = n_per_sess;
  !!   nc_mrds = n_traps_per_sess;
  !! } else {
  !!   nr_mrds = 1;
  !!   nc_mrds = 1;
  !! }
  // Expected supplementary information from each mask point.
  init_3darray mrds_dist(1,n_sessions,1,nr_mrds,1,nc_mrds)
  init_3darray dists(1,n_sessions,1,n_traps_per_sess,1,n_mask_per_sess)
  init_3darray angs(1,n_sessions,1,nr_angmat,1,nc_angmat)
  init_3darray toa_ssq(1,n_sessions,1,nr_toa_ssq,1,nc_toa_ssq)
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
  // Sorting out density covariates.
  ivector nr_ihd(1,n_sessions)
  int nc_ihd
  int length_local_D_mask
  !!   nr_ihd = n_mask_per_sess;
  !!   nc_ihd = n_D_betapars;
  init_3darray mm_ihd(1,n_sessions,1,nr_ihd,1,nc_ihd)
  // Setting n_local permanently if local integration is disabled.
  int nr_localmats
  int nc_localmats
  int n_local
  number dir
  number bearing_to_trap
  number orientation

PARAMETER_SECTION
  objective_function_value f
  // Some sort of individual contribution to the likelihood.
  number f_ind
  // Something-or-other required for directional calling.
  matrix fs(1,n_sessions,1,length_fs)
  // Linked density parameters.
  init_bounded_number_vector D_betapars_link(1,n_D_betapars,D_betapars_lb,D_betapars_ub,D_betapars_phase)
  // Linked detection function parameters.
  init_bounded_number_vector detpars_link(1,n_detpars,detpars_lb,detpars_ub,detpars_phase)
  // Linked supplementary information parameters.
  init_bounded_number_vector suppars_link(1,n_suppars,suppars_lb,suppars_ub,suppars_phase)
  vector D_ihdpars(1,n_D_betapars+1);
  !! D_betapars_link.set_scalefactor(D_betapars_sf);
  !! detpars_link.set_scalefactor(detpars_sf);
  !! suppars_link.set_scalefactor(suppars_sf);
  // I don't know why this scalefactor has to be so big but it just seems to work OK?
  //!! D_betapars.set_scalefactor(10000000); Uh I'm commenting this out because the comment above is weird.
  // Collection of parameter estimates.
  sdreport_vector par_ests(1,n_ests)
  // Effective sampling area.
  sdreport_vector esa(1,n_sessions)
  // Log-probability of capture of animals at each mask point by each detector.
  4darray log_capt_probs(1,n_sessions,1,dummy_n_dir_quadpoints,1,n_traps_per_sess,1,n_mask_per_sess)
  // Log-probability of evasion of animals at each mask point by each detector.
  4darray log_evade_probs(1,n_sessions,1,dummy_n_dir_quadpoints,1,n_traps_per_sess,1,n_mask_per_sess)
  // Expected received signal strength at each detector by a sound emitted at each mask point.
  4darray expected_ss(1,n_sessions,1,dummy_n_dir_quadpoints,1,nr_expected_ss,1,nc_expected_ss)
  // Some parameters or something.
  sdreport_number D
  number corr_ss
  number cond_corr_ss
  vector detpars(1,n_detpars)
  vector suppars(1,n_suppars)
  // It seems like you don't need the following line, but if you
  // delete it you get a compilation error on Macs. NOTE: I've
  // commented it out because with multi-session models n_traps is not
  // constant. If compilation fails on Mac this might be why. Filippo
  // reported no problems so it's probably fine!
  //matrix capt_hist(1,n_sessions,1,n_traps)
  // A bunch of variables that are used for stuff.
  number det_prob
  vector sum_det_probs(1,n_sessions)
  vector sum_D_det_probs(1,n_sessions)
  matrix D_mask(1,n_sessions,1,nr_ihd)
  vector D_mask_vec(1,n_mask_total)
  number undet_prob
  number undet_lower_prob
  number capt_prob
  number ss_resid
  number point_capt
  number point_evade
  number diag_sigma_ss
  number offdiag_sigma_ss

PROCEDURE_SECTION
  // Grabbing detection function.
  detfn_pointer detfn = get_detfn(detfn_id);
  // Converting linked parameters to real parameters and setting up par_ests vector.
  j = 1;
  for (i = 1; i <= n_D_betapars; i++){
    par_ests(j) = D_betapars_link(i);
    j++;
  }
  // Getting D for homogeneous density models.
  D = mfexp(D_betapars_link(1));
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
  // Generating variance-covariance matrix and correlation matrix for fits with heterogeneity in source strength.
  if (fit_het_source){
    diag_sigma_ss = square(detpars(4)) + square(detpars(5));
    offdiag_sigma_ss = square(detpars(4));
    corr_ss = square(detpars(4))/(square(detpars(4)) + square(detpars(5)));
  }
  // Start of likelihood calculation.
  f = 0.0;
  fs = 0.0;
  // Calculating density at each mask point.
  for (s = 1; s <= n_sessions; s++){
    D_mask.rowfill(s, mfexp(mm_ihd(s)*D_betapars_link));
  }  
  i = 1;
  for (s = 1; s <= n_sessions; s++){
    for (j = 1; j <= n_mask_per_sess(s); j++){
      D_mask_vec(i) = D_mask(s, j);
      i++;
    }
  } 
  // Calculating mask detection probabilities and expected signal strengths...
  sum_det_probs = 0;
  sum_D_det_probs = 0;
  // ... for a directional model...
  if (fit_dir){
    for (s = 1; s <= n_sessions; s++){
      for (b = 1; b <= n_dir_quadpoints; b++){
        dir = 2*M_PI*(b - 1)/n_dir_quadpoints;
        for (i = 1; i <= n_mask_per_sess(s); i++){
          point_evade = 1;
          for (j = 1; j <= n_traps_per_sess(s); j++){
            bearing_to_trap = angs(s, j, i);
            // Adjusting bearing_to_trap, angs matrix gives bearing from
            // trap to animal. Guess this doesn't really matter if
            // f(b_i) is uniform on [0, 2*pi), but nice to have correct.
            if (bearing_to_trap > M_PI){
              bearing_to_trap -= M_PI;
            } else {
              bearing_to_trap += M_PI;
            }
            // Orientation of trap with respect to direction of animal.
            orientation = bearing_to_trap - dir;
            dist = dists(s, j, i);
            capt_prob = detfn(dist, detpars, cutoff, orientation);
            log_capt_probs(s, b, j, i) = log(capt_prob + dbl_min);
            log_evade_probs(s, b, j, i) = log(1 - capt_prob + dbl_min);
            point_evade *= 1 - capt_prob;
            if (fit_ss){
              if (linkfn_id == 3){
                expected_ss(s, b, j, i) = detpars(1) - 10*log10(square(dist)) - (detpars(2) - (detpars(3)*(cos(orientation) - 1.0)/2))*(dist - 1.0);
              } else {
                expected_ss(s, b, j, i) = detpars(1) - (detpars(2) - (detpars(3)*(cos(orientation) - 1.0)/2))*dist;	      
                if (linkfn_id == 2){
                  expected_ss(s, b, j, i) = mfexp(expected_ss(s, b, j, i));
                }
              }
            }
          }
          point_capt = 1 - point_evade;
          sum_det_probs(s) += point_capt/n_dir_quadpoints;
	  sum_D_det_probs(s) += D_mask(s, i)*point_capt/n_dir_quadpoints;
        }
      }
    }
  // ... and a non-directional model.
  } else {
    orientation = 0;
    for (s = 1; s <= n_sessions; s++){
      for (i = 1; i <= n_mask_per_sess(s); i++){
        // Saving expected signal strengths.
        if (fit_ss){
          dvar_vector mu_ss(1, n_traps_per_sess(s));
          if (linkfn_id == 3){
            mu_ss = detpars(1) - 10*log10(square(column(dists(s), i))) - detpars(2)*(column(dists(s), i) - 1.0);
          } else {
            mu_ss = detpars(1) - detpars(2)*column(dists(s), i);
            if (linkfn_id == 2){
              mu_ss = mfexp(mu_ss);
            }
          }  
          expected_ss(s, 1).colfill(i, mu_ss);
          // For a model with heterogeneity in signal strengths.
          if (fit_het_source){
            dvar_vector z_ss(1, n_traps_per_sess(s));
            z_ss = (cutoff - mu_ss)/pow(diag_sigma_ss, 0.5);
            undet_prob = pmvn(z_ss, corr_ss, het_source_gh, het_source_weights, het_source_nodes, n_het_source_quadpoints, -5, 5);
          }
        }
        if (!fit_het_source){
          // No heterogeneity in source signal strengths.
          undet_prob = 1;
          for (j = 1; j <= n_traps_per_sess(s); j++){
            dist = dists(s, j, i);
            capt_prob = detfn(dist, detpars, cutoff, orientation);
            // Compare to calculating these outside loop.
            log_capt_probs(s, 1, j, i) = log(capt_prob + dbl_min);
            log_evade_probs(s, 1, j, i) = log(1 - capt_prob + dbl_min);
            undet_prob *= 1 - capt_prob;
          }
        }
        det_prob = 1 - undet_prob;
        sum_det_probs(s) += det_prob + dbl_min;
        sum_D_det_probs(s) += D_mask(s, i)*det_prob + dbl_min;
      }
    }
  }
  // Contribution due to capture history.
  for (s = 1; s <= n_sessions; s++){
    for (u = 1; u <= n_unique_per_sess(s); u++){
      dvar_vector capt_hist(1,n_traps_per_sess(s));
      capt_hist = row(capt_bin_unique(s), u);
      if (local == 1){
        n_local = n_local_per_unique(s, u);
	nr_localmats = n_traps_per_sess(s);
        nc_localmats = n_local;
      } else {
      	n_local = n_mask_per_sess(s);
	nr_localmats = 1;
	nc_localmats = 1;
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
            local_angs.colfill(j, column(angs(s), which_local_per_unique(s, u, j)));
          }
          angs_pointer = &local_angs;
        } else {
          angs_pointer = &angs(s);
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
            local_dists.colfill(j, column(dists(s), which_local_per_unique(s, u, j)));
          }
          dists_pointer = &local_dists;
        } else {
          dists_pointer = &dists(s);
        }
      }
      // Filling local_toa_ssq.
      if (fit_toas){
        nr_local_toa_ssq = n_per_sess(s);
        nc_local_toa_ssq = nc_localmats;
      } else {
        nr_local_toa_ssq = 1;
        nc_local_toa_ssq = 1;
      }
      dmatrix * toa_ssq_pointer;
      dmatrix local_toa_ssq(1, nr_local_toa_ssq(s), 1, nc_local_toa_ssq(s));
      if (fit_toas){
        if (local){
          for (j = 1; j <= n_local; j++){
            local_toa_ssq.colfill(j, column(toa_ssq(s), which_local_per_unique(s, u, j)));
          }
          toa_ssq_pointer = &local_toa_ssq;
        } else {
          toa_ssq_pointer = &toa_ssq(s);
        }
      }
      // If dealing with inhomogeneous density, filling local densities.
      length_local_D_mask = n_local;
      dvar_vector * D_mask_pointer;
      dvar_vector local_D_mask(1, length_local_D_mask);
      dvar_vector D_contrib(1, length_local_D_mask);
      if (local){
        for (j = 1; j <= n_local; j++){
          local_D_mask(j) = D_mask(s, which_local_per_unique(s, u, j));
        }
        D_mask_pointer = &local_D_mask;
      } else {
        D_mask_pointer = &D_mask(s);
      }
      D_contrib = log(*D_mask_pointer + dbl_min);
      // Getting capture probabilities at local mask points.
      dvar_matrix * log_capt_probs_pointer;
      dvar_matrix * log_evade_probs_pointer;
      dvar_matrix local_log_capt_probs(1,nr_localmats,1,nc_localmats);
      dvar_matrix local_log_evade_probs(1,nr_localmats,1,nc_localmats);
      for (b = 1; b <= n_dir_quadpoints; b++){
        // Resetting i index.
        i = 0;
        if (local){
          for (j = 1; j <= n_local; j++){
            local_log_capt_probs.colfill(j, column(log_capt_probs(s, b), which_local_per_unique(s, u, j)));
            local_log_evade_probs.colfill(j, column(log_evade_probs(s, b), which_local_per_unique(s, u, j)));
          }
          log_capt_probs_pointer = &local_log_capt_probs;
          log_evade_probs_pointer = &local_log_evade_probs;
        } else {
          log_capt_probs_pointer = &log_capt_probs(s, b);
          log_evade_probs_pointer = &log_evade_probs(s, b);
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
        dvar_matrix local_expected_ss(1,nr_local_expected_ss(s),1,nc_local_expected_ss(s));
        if (fit_ss){
          if (local){
            for (j = 1; j <= n_local; j++){
              local_expected_ss.colfill(j, column(expected_ss(s, b), which_local_per_unique(s, u, j)));
            }
            expected_ss_pointer = &local_expected_ss;
          } else {
            expected_ss_pointer = &expected_ss(s, b);
          }
        }
        i_start = 0;
        k = 1;
        while (k < u){
          i_start += capt_bin_freqs(s, k);
          k++;
        }
        i_start += 1;
        i_end = i_start + capt_bin_freqs(s, u) - 1;
        for (i = i_start; i <= i_end; i++){
          // Contribution from capture locations.
          if (fit_ss){
            // Do something in here for heterogeneous source strengths.
            if (fit_het_source){
              double n_dets = value(sum(capt_hist));
              bool all_dets = n_dets == n_traps_per_sess(s);
              // Observed signal strengths.
              dvector obs_ss(1, n_dets);
              if (all_dets){
                obs_ss = row(capt_ss(s), i);
                for (j = 1; j <= n_local; j++){
                  bincapt_contrib(j) = log_dmvn_diag(obs_ss, column((*expected_ss_pointer), j), diag_sigma_ss, offdiag_sigma_ss, dbl_min);
                }
                // Put in log_dmvn for received signal strengths.
              } else {
                double n_nodets = n_traps_per_sess(s) - n_dets;
                ivector ind_det(1, n_dets);
                ivector ind_nodet(1, n_nodets);
                a = 1;
                c = 1;
                for (j = 1; j <= n_traps_per_sess(s); j++){
                  if (capt_hist(j) == 1){
                    ind_det(a) = j;
                    a++;
                  } else {
                    ind_nodet(c) = j;
                    c++;
                  }
                }
                obs_ss = row(capt_ss(s), i)(ind_det);
                // Marginal and conditional mean vectors.
	        dvar_vector mu_ss_det(1, n_dets);
                dvar_vector mu_ss_nodet(1, n_nodets);
                dvar_vector mu_ss_cond(1, n_nodets);
                // Diagonal and off-diagonal elements of the conditional variance-covariance matrix.
                dvariable diag_sigma_ss_cond = (pow(detpars(5), 4) + (n_dets + 1)*square(detpars(5))*square(detpars(4)))/(square(detpars(5)) + n_dets*square(detpars(4)));
                dvariable offdiag_sigma_ss_cond = (square(detpars(5))*square(detpars(4)))/(square(detpars(5)) + n_dets*square(detpars(4)));
                // Conditional correlation.
                dvariable corr_ss_cond = square(detpars(4))/(square(detpars(5)) + (n_dets + 1)*square(detpars(4)));
                // Standardised discrepency between conditional expected strength and cutoff strength.
                dvar_vector z_ss_nodet(1, n_nodets);
                // Joint density of capture history and something else (source strength?), conditional on location.
                for (j = 1; j <= n_local; j++){
                  mu_ss_det = column((*expected_ss_pointer), j)(ind_det);
                  mu_ss_nodet = column((*expected_ss_pointer), j)(ind_nodet);
                  mu_ss_cond = mu_ss_nodet + (square(detpars(4))*sum(obs_ss - mu_ss_det))/(square(detpars(5)) + n_dets*square(detpars(4)));
                  bincapt_contrib(j) = log_dmvn_diag(obs_ss, mu_ss_det, diag_sigma_ss, offdiag_sigma_ss, dbl_min);
                  z_ss_nodet = (cutoff - mu_ss_cond)/pow(diag_sigma_ss_cond, 0.5);
                  bincapt_contrib(j) += log(pmvn(z_ss_nodet, corr_ss_cond, het_source_gh, het_source_weights, het_source_nodes, n_het_source_quadpoints, -5, 5) + dbl_min);
                }
              }
            } else {
              // I think this part could be improved by only calculating
	      // local_log_ss_density for detected traps?
              dvar_matrix local_log_ss_density(1,n_traps_per_sess(s),1,n_local);
              for (j = 1; j <= n_traps_per_sess(s); j++){
                local_log_ss_density.rowfill(j, log_dnorm(capt_ss(s, i, j), row((*expected_ss_pointer), j), detpars(5)));
              }
              bincapt_contrib = capt_hist*local_log_ss_density + evade_contrib;
            }
          }
          if (any_suppars){
            dvar_vector supp_contrib(1,n_local);
            supp_contrib = 0;
            for (j = 1; j <= n_traps_per_sess(s); j++){
              //  Try setting up a ragged array of capture locations for each individual instead.
              if (capt_bin_unique(s, u, j)){
                // Contribution from bearings.
                if (fit_angs){
                  supp_contrib += suppars(kappa_ind)*cos(capt_ang(s, i, j) - row((*angs_pointer), j));
                }
                // Contribution from distances.
                if (fit_dists){  
                  // Try saving alpha separately.
	          dvar_vector beta(1, n_local);
	          beta = suppars(alpha_ind)/row((*dists_pointer), j);
                  supp_contrib += suppars(alpha_ind)*log(beta) + (suppars(alpha_ind) - 1)*log(capt_dist(s, i, j)) - beta*capt_dist(s, i, j);
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
	    f_ind = sum(mfexp(bincapt_contrib + supp_contrib + D_contrib));
          } else {
            f_ind = sum(mfexp(bincapt_contrib + D_contrib));
          }
          // For directional calling, save components due to each
          // direction for each individual.
          if (fit_dir){
            fs(s, i) += f_ind/n_dir_quadpoints;
          } else {
            f -= log(f_ind + dbl_min);
          }
	  //cout << "f_ind(" << i << "): " << log(f_ind) << endl;
        }
      }
    }
  }
  // For directional calling, fs contains liklihood contributions for
  // each individual. Need to log and sum for log-likelihood.
  if (fit_dir){
    f -= sum(log(fs + dbl_min));
  }
  for (s = 1; s <= n_sessions; s++){
    // Calculating ESAs.
    esa(s) = A_per_sess(s)*sum_det_probs(s);
    // Adding contribution from ns.
    f -= log_dpois(n_per_sess(s), survey_length(s)*A_per_sess(s)*sum_D_det_probs(s));
    // Extra bit that falls out of log-likelihood.
    f -= -n_per_sess(s)*log(sum_D_det_probs(s));
  }
  // Printing trace.
  if (trace){
    for (i = 1; i <= n_D_betapars; i++){
      cerr << "D par " << i << ": " << D_betapars_link(i) << ", ";
    }
    for (i = 1; i <= n_detpars; i++){
      cerr << "DF Par " << i << ": " << detpars(i) << ", ";
    }
    if (any_suppars){
      for (i = 1; i <= n_suppars; i++){
        cerr << "Supp Par " << i << ": " << suppars(i) << ", ";
      }
    }
    cerr << "LL: " << -f << ", ESA: " << esa << endl;
  }

REPORT_SECTION
  // Writing D to report file.
  report << "# D:" << endl << D << endl;
  // Writing ESAs to report file.
  for (i = 1; i <= n_sessions; i++){
    report << "# esa[" << i << "]:" << endl << esa(i) << endl;
  }
  // Writing D_mask to report file.
  for (i= 1; i <= n_sessions; i++){
    report << "# D_mask[" << i << "]:" << endl << D_mask(i) << endl;
  }
  
GLOBALS_SECTION
  #include <detfuns.cpp>
  #include <helpers.cpp>
  #include <invlinkfuns.cpp>
  #include <densfuns.cpp>

