#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

double log_dnorm (double x, double mu, const double sigma)
{
  return -0.5*log(2*M_PI) - log(sigma) - pow(x - mu, 2)/(2*pow(sigma, 2));
}


//' Evaluating the likelihood using C++
//'
//' Returns the likelihood for a vector of parameter values.
//'
//' @param pars Parameter values.
//' @param objs Various other objects.
//'
//' @return The value of the negative log-likelihood.
//'
//' @export
// [[Rcpp::export]]
double secr_nll(const NumericVector& link_pars, const List& dat){
  double D = exp(link_pars[0]);
  double b0_ss = exp(link_pars[1]);
  double b1_ss = exp(link_pars[2]);
  double sigma_ss = exp(link_pars[3]);
  int n_unique = as<int>(dat["n_unique"]);
  int n = as<int>(dat["n"]);
  int n_traps = as<int>(dat["n_traps"]);
  int n_mask = as<int>(dat["n_mask"]);
  IntegerMatrix capt_bin_unique = as<IntegerMatrix>(dat["capt_bin_unique"]);
  IntegerVector capt_bin_freqs = as<IntegerVector>(dat["capt_bin_freqs"]);
  int first_calls = as<int>(dat["first_calls"]);
  double cutoff = as<double>(dat["cutoff"]);
  NumericVector cutoff_vec(1);
  cutoff_vec[0] = cutoff;
  double lower_cutoff = as<double>(dat["lower_cutoff"]);
  NumericVector lower_cutoff_vec(1);
  lower_cutoff_vec[0] = lower_cutoff;
  NumericMatrix capt_ss = as<NumericMatrix>(dat["capt_ss"]);
  NumericMatrix dists = as<NumericMatrix>(dat["dists"]);
  int trace = as<int>(dat["trace"]);
  double A = as<double>(dat["A"]);
  double dbl_min = as<double>(dat["DBL_MIN"]);
  int i, j;
  double mu_ss;
  double det_prob;
  double sub_det_prob;
  double undet_prob;
  double undet_lower_prob;
  double sum_det_probs = 0;
  double sum_sub_det_probs = 0;
  NumericVector mask_det_probs(n_mask);
  NumericVector mask_all_det_probs(n_mask);
  NumericMatrix expected_ss(n_traps, n_mask);
  NumericMatrix log_capt_probs(n_traps, n_mask);
  NumericMatrix log_evade_probs(n_traps, n_mask);
  double capt_prob;
  for (i = 0; i < n_mask; i++){
    undet_prob = 1;
    undet_lower_prob = 1;
    for (j = 0; j < n_traps; j++){
      mu_ss = b0_ss - b1_ss*dists(j, i);
      expected_ss(j, i) = mu_ss;
      // For some reason first argument to pnorm must be a vector.
      capt_prob = 1 - pnorm(cutoff_vec, mu_ss, sigma_ss)[0];
      log_capt_probs(j, i) = log(capt_prob + dbl_min);
      log_evade_probs(j, i) = log(1 - capt_prob + dbl_min);
      undet_prob *= 1 - capt_prob;
      if (first_calls){
	// As above.
        undet_lower_prob *= pnorm(lower_cutoff_vec, mu_ss, sigma_ss)[0];
      }
    }
    det_prob = 1 - undet_prob;
    sum_det_probs += det_prob;
    if (first_calls){
      // Saving detection probabilities for first calls, and
      // detection probabilities for any subsequent calls.
      mask_det_probs[i] = 0;
      mask_all_det_probs[i] = 0;
      if (1 - undet_lower_prob > 1e-30){
	sum_sub_det_probs += (det_prob*undet_lower_prob)/(1 - undet_lower_prob);
	sub_det_prob = (det_prob*undet_lower_prob)/(1 - undet_lower_prob);
	mask_det_probs[i] += det_prob;
	mask_all_det_probs[i] += det_prob + sub_det_prob;
      }
    }
  }
  NumericVector capt_hist(n_traps);
  NumericVector evade_contrib(n_mask);
  NumericVector bincapt_contrib(n_mask);
  NumericMatrix log_ss_density(n_traps, n_mask);
  double f = 0.0;
  double f_ind;
  int i_start, i_end, m, u, k;
  // Contribution from capture history.
  for (u = 0; u < n_unique; u++){
    capt_hist = capt_bin_unique(u, _);
    // Calculating evasion contribution.
    for (m = 0; m < n_mask; m++){
      evade_contrib[m] = 0;
      for (j = 0; j < n_traps; j++){
	if (capt_hist[j] == 0){
	  evade_contrib[m] += log_evade_probs(j, m);
	}
      }
    }
    // Sorting indices for individuals with uth unique capture history.
    i_start = -1;
    k = 0;
    while (k < u){
      i_start += capt_bin_freqs[k];
      k++;
    }
    i_start += 1;
    i_end = i_start + capt_bin_freqs[u];
    for (i = i_start; i < i_end; i++){
      for (m = 0; m < n_mask; m++){
	bincapt_contrib[m] = evade_contrib[m];
	for (j = 0; j < n_traps; j++){
	  if (capt_hist[j] == 1){
	    bincapt_contrib[m] += log_dnorm(capt_ss(i, j), expected_ss(j, m), sigma_ss);
	  }
	}
      }
      f_ind = sum(exp(bincapt_contrib + log(mask_all_det_probs + dbl_min) - log(mask_det_probs + dbl_min)));
      f -= log(f_ind + dbl_min);
    }
  }
  double esa = A*(sum_det_probs + sum_sub_det_probs);
  NumericVector n_vec(1);
  n_vec[0] = n;
  f -= log(dpois(n_vec, D*esa)[0]);
  f -= -n*log(sum_det_probs + sum_sub_det_probs);
  if (trace){
    cout << "D: " << D << ", b0.ss: " << b0_ss << ", b1.ss: " << b1_ss << ", sigma.ss: " << sigma_ss << ", LL: " << -f << endl;
  }
  return f;
}

