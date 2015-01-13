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
  int b
  int a
  int c
  int d
  int i_start
  int i_end
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
  // Logical indicator for directional calling.
  init_int fit_dir
  init_number n_dir_quadpoints
  init_int fit_het_source
  init_int het_source_gh
  init_number n_het_source_quadpoints
  int n_het_source_nodes
  !! if (fit_het_source & het_source_gh){
  !!   n_het_source_nodes = n_het_source_quadpoints;
  !! } else {
  !!   n_het_source_nodes = 1;
  !! }
  init_vector het_source_nodes(1,n_het_source_nodes);
  init_vector het_source_weights(1,n_het_source_nodes);
  int length_fs
  !! if (fit_dir){
  !!   length_fs = n;
  !! } else {
  !!   length_fs = 1;
  !! }
  int nr_ang
  int nc_ang
  int nr_angmat
  int nc_angmat
  int nr_local_angmat
  int nc_local_angmat
  !! if (fit_angs){
  !!   nr_ang = n;
  !!   nc_ang = n_traps;
  !!   nr_angmat = n_traps;
  !!   nc_angmat = n_mask;
  !! } else if (fit_dir) {
  !!   nr_ang = 1;
  !!   nc_ang = 1;
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
  init_int first_calls
  int n_mask_det_probs
  !! if (first_calls == 1){
  !!   n_mask_det_probs = n_mask;
  !! } else {
  !!   n_mask_det_probs = 1;
  !! }
  init_number lower_cutoff
  init_number first_calls_trunc
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
  number dir
  number bearing_to_trap
  number orientation
  vector quo(1,n_mask)
  vector den(1,n_mask)

PARAMETER_SECTION
  objective_function_value f
  number f_ind
  vector fs(1,length_fs)
  init_bounded_number D_link(D_lb,D_ub,D_phase)
  init_bounded_number_vector detpars_link(1,n_detpars,detpars_lb,detpars_ub,detpars_phase)
  init_bounded_number_vector suppars_link(1,n_suppars,suppars_lb,suppars_ub,suppars_phase)
  !! D_link.set_scalefactor(D_sf);
  !! detpars_link.set_scalefactor(detpars_sf);
  !! suppars_link.set_scalefactor(suppars_sf);
  sdreport_vector par_ests(1,n_ests)
  sdreport_number esa
  3darray log_capt_probs(1,n_dir_quadpoints,1,n_traps,1,n_mask)
  3darray log_evade_probs(1,n_dir_quadpoints,1,n_traps,1,n_mask)
  3darray expected_ss(1,n_dir_quadpoints,1,nr_expected_ss,1,nc_expected_ss)
  number D
  number corr_ss
  number cond_corr_ss
  vector detpars(1,n_detpars)
  vector suppars(1,n_suppars)
  // Keep this here, otherwise get an error on Mac.
  vector capt_hist(1,n_traps)
  // Need to save probabilities of detection and probabilities of subsequent detection for first call models.
  vector mask_det_probs(1,n_mask_det_probs)
  vector mask_all_det_probs(1,n_mask_det_probs)
  number det_prob
  number sum_det_probs
  // Probabilities for a subsequent call appearing in capture history.
  number sub_det_prob
  number sum_sub_det_probs
  number undet_prob
  // Probabilities of a call being undetected at the lower cutoff.
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
  // Generating variance-covariance matrix and correlation matrix for fits with heterogeneity in source strength.
  if (fit_het_source){
    diag_sigma_ss = square(detpars(4)) + square(detpars(5));
    offdiag_sigma_ss = square(detpars(4));
    corr_ss = square(detpars(4))/(square(detpars(4)) + square(detpars(5)));
  }
  // Start of likelihood calculation.
  f = 0.0;
  fs = 0.0;
  // Setting up sum_sub_det_probs for a first call only model.
  sum_sub_det_probs = 0;
  // Calculating mask detection probabilities and expected signal strengths...
  sum_det_probs = 0;
  // ... for a directional model...
  if (fit_dir){
    for (b = 1; b <= n_dir_quadpoints; b++){
      dir = 2*M_PI*(b - 1)/n_dir_quadpoints;
      for (i = 1; i <= n_mask; i++){
        point_evade = 1;
        for (j = 1; j <= n_traps; j++){
          bearing_to_trap = angs(j, i);
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
          dist = dists(j, i);
          capt_prob = detfn(dist, detpars, cutoff, orientation);
          log_capt_probs(b, j, i) = log(capt_prob + DBL_MIN);
          log_evade_probs(b, j, i) = log(1 - capt_prob + DBL_MIN);
          point_evade *= 1 - capt_prob;
          if (fit_ss){
            if (linkfn_id == 3){
              expected_ss(b, j, i) = detpars(1) - 10*log10(square(dist)) - (detpars(2) - (detpars(3)*(cos(orientation) - 1.0)))*(dist - 1.0);
            } else {
              expected_ss(b, j, i) = detpars(1) - (detpars(2) - (detpars(3)*(cos(orientation) - 1)))*dist;
              if (linkfn_id == 2){
                expected_ss(b, j, i) = mfexp(expected_ss(b, j, i));
              }
            }
          }
        }
        point_capt = 1 - point_evade;
        sum_det_probs += point_capt/n_dir_quadpoints;
      }
    }
  // ... and a non-directional model.
  } else {
    orientation = 0;
    for (i = 1; i <= n_mask; i++){
      // Saving expected signal strengths.
      if (fit_ss){
        dvar_vector mu_ss(1, n_traps);
        if (linkfn_id == 3){
          mu_ss = detpars(1) - 10*log10(square(column(dists, i))) - detpars(2)*(column(dists, i) - 1.0);
        } else {
          mu_ss = detpars(1) - detpars(2)*column(dists, i);
          if (linkfn_id == 2){
            mu_ss = mfexp(mu_ss);
          }
        }  
        expected_ss(1).colfill(i, mu_ss);
        // For a model with heterogeneity in signal strengths.
        if (fit_het_source){
          dvar_vector z_ss(1, n_traps);
          z_ss = (cutoff - mu_ss)/pow(diag_sigma_ss, 0.5);
          undet_prob = pmvn(z_ss, corr_ss, het_source_gh, het_source_weights, het_source_nodes, n_het_source_quadpoints, -5, 5);
        }
      }
      if (!fit_het_source){
        // No heterogeneity in source signal strengths.
        undet_prob = 1;
	undet_lower_prob = 1;
        for (j = 1; j <= n_traps; j++){
          dist = dists(j, i);
          capt_prob = detfn(dist, detpars, cutoff, orientation);
          // Compare to calculating these outside loop.
          log_capt_probs(1, j, i) = log(capt_prob + DBL_MIN);
          log_evade_probs(1, j, i) = log(1 - capt_prob + DBL_MIN);
          undet_prob *= 1 - capt_prob;
	  if (first_calls){
            undet_lower_prob *= 1 - detfn(dist, detpars, lower_cutoff, orientation);
          }
        }
      }
      det_prob = 1 - undet_prob;
      sum_det_probs += det_prob + DBL_MIN;
      if (first_calls){
	// Saving detection probabilities for first calls, and
	// detection probabilities for any subsequent calls.
	mask_det_probs(i) = 0;
	mask_all_det_probs(i) = 0;
        //num = det_prob*undet_lower_prob;
	quo(i) = value((det_prob*undet_lower_prob)/(1 - undet_lower_prob + DBL_MIN));
	den(i) = value(1 - undet_lower_prob);
	//if (den(i) > first_calls_trunc){
	if (den(i) > 1e-20){
          // Will probably find some numerical instability here.
          sum_sub_det_probs += (det_prob*undet_lower_prob)/(1 - undet_lower_prob);
	}
	if (den(i) > 1e-20){
          sub_det_prob = (det_prob*undet_lower_prob)/(1 - undet_lower_prob);
	  mask_det_probs(i) += det_prob;
	  mask_all_det_probs(i) += det_prob + sub_det_prob;
        }
      }
    }  
  }
  // Contribution due to capture history.
  for (u = 1; u <= n_unique; u++){
    capt_hist = row(capt_bin_unique, u);
    if (local == 1){
      n_local = all_n_local(u);
      nc_localmats = n_local;
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
          local_log_capt_probs.colfill(j, column(log_capt_probs(b), all_which_local(u, j)));
          local_log_evade_probs.colfill(j, column(log_evade_probs(b), all_which_local(u, j)));
        }
        log_capt_probs_pointer = &local_log_capt_probs;
        log_evade_probs_pointer = &local_log_evade_probs;
      } else {
        log_capt_probs_pointer = &log_capt_probs(b);
        log_evade_probs_pointer = &log_evade_probs(b);
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
            local_expected_ss.colfill(j, column(expected_ss(b), all_which_local(u, j)));
          }
          expected_ss_pointer = &local_expected_ss;
        } else {
          expected_ss_pointer = &expected_ss(b);
        }
      }
      i_start = 0;
      k = 1;
      while (k < u){
        i_start += capt_bin_freqs(k);
        k++;
      }
      i_start += 1;
      i_end = i_start + capt_bin_freqs(u) - 1;
      for (i = i_start; i <= i_end; i++){
        // Contribution from capture locations.
        if (fit_ss){
          // Do something in here for heterogeneous source strengths.
          if (fit_het_source){
            double n_dets = value(sum(capt_hist));
            bool all_dets = n_dets == n_traps;
            // Observed signal strengths.
            dvector obs_ss(1, n_dets);
            if (all_dets){
              obs_ss = row(capt_ss, i);
              for (j = 1; j <= n_local; j++){
                bincapt_contrib(j) = log_dmvn_diag(obs_ss, column((*expected_ss_pointer), j), diag_sigma_ss, offdiag_sigma_ss, DBL_MIN);
              }
              // Put in log_dmvn for received signal strengths.
            } else {
              double n_nodets = n_traps - n_dets;
              ivector ind_det(1, n_dets);
              ivector ind_nodet(1, n_nodets);
              a = 1;
              c = 1;
              for (j = 1; j <= n_traps; j++){
                if (capt_hist(j) == 1){
                  ind_det(a) = j;
                  a++;
                } else {
                  ind_nodet(c) = j;
                  c++;
                }
              }
              obs_ss = row(capt_ss, i)(ind_det);
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
              // In here goes f(w, r | x), right?
              //
              // Why are you writing questions? I don't fucking know
              // but it seems to work so I guess you did this fine.
              for (j = 1; j <= n_local; j++){
                mu_ss_det = column((*expected_ss_pointer), j)(ind_det);
                mu_ss_nodet = column((*expected_ss_pointer), j)(ind_nodet);
                mu_ss_cond = mu_ss_nodet + (square(detpars(4))*sum(obs_ss - mu_ss_det))/(square(detpars(5)) + n_dets*square(detpars(4)));
                bincapt_contrib(j) = log_dmvn_diag(obs_ss, mu_ss_det, diag_sigma_ss, offdiag_sigma_ss, DBL_MIN);
                z_ss_nodet = (cutoff - mu_ss_cond)/pow(diag_sigma_ss_cond, 0.5);
                bincapt_contrib(j) += log(pmvn(z_ss_nodet, corr_ss_cond, het_source_gh, het_source_weights, het_source_nodes, n_het_source_quadpoints, -5, 5) + DBL_MIN);
              }
            }
          } else {
	    // I think this part could be improved by only calculating
	    // local_log_ss_density for detected traps?
            dvar_matrix local_log_ss_density(1,n_traps,1,n_local);
            for (j = 1; j <= n_traps; j++){
              local_log_ss_density.rowfill(j, log_dnorm(capt_ss(i, j), row((*expected_ss_pointer), j), detpars(5)));
            }
            bincapt_contrib = capt_hist*local_log_ss_density + evade_contrib;
          }
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
	  if (first_calls){
            //f_ind = sum(mfexp(bincapt_contrib + supp_contrib + log_mask_all_det_probs - log_mask_det_probs));
	    f_ind = sum(mfexp(bincapt_contrib + supp_contrib + log(mask_all_det_probs + DBL_MIN) - log(mask_det_probs + DBL_MIN)));
          } else {
            f_ind = sum(mfexp(bincapt_contrib + supp_contrib));
          }
        } else {
	  if (first_calls){
            //f_ind = sum(mfexp(bincapt_contrib + log_mask_all_det_probs - log_mask_det_probs));
	    f_ind = sum(mfexp(bincapt_contrib + log(mask_all_det_probs + DBL_MIN) - log(mask_det_probs + DBL_MIN)));
          } else {
            f_ind = sum(mfexp(bincapt_contrib));
          }
        }
        // For directional calling, save components due to each
        // direction for each individual.
        if (fit_dir){
          fs(i) += f_ind/n_dir_quadpoints;
        } else {
          f -= log(f_ind + DBL_MIN);
        }
      }
    }
  }
  // For directional calling, fs contains liklihood contributions for
  // each individual. Need to log and sum for log-likelihood.
  if (fit_dir){
    f -= sum(log(fs + DBL_MIN));
  }
  // Calculating ESA.
  esa = A*(sum_det_probs + sum_sub_det_probs);
  // Contribution from n.
  f -= log_dpois(n, D*esa);
  // Extra bit that falls out of ll.
  f -= -n*log(sum_det_probs + sum_sub_det_probs);
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
    cout << "LL: " << -f << ", ESA: " << esa << ", test" << endl;
  }

REPORT_SECTION
  // Writing ESA to report file.
  report << esa << endl;
  // Writing subsequent calls component of ESA to report file.
  report << sum_sub_det_probs << endl;
  // Writing quotients and denominators for first calls to report file.
  report << den << endl;
  report << quo << endl;

GLOBALS_SECTION
  #include <detfuns.cpp>
  #include <helpers.cpp>
  #include <invlinkfuns.cpp>
  #include <densfuns.cpp>

