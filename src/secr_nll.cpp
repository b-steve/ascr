#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

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
double secr_nll(const NumericVector& pars, const List& dat){
  double D = exp(pars[0]);
  double b0_ss = exp(pars[1]);
  double b1_ss = exp(pars[2]);
  double sigma_ss = exp(pars[3]);
  int n_unique = as<int>(dat["n_unique"]);
  int n = as<int>(dat["n"]);
  int n_traps = as<int>(dat["n_traps"]);
  int n_mask = as<int>(dat["n_mask"]);
  IntegerMatrix capt_bin_unique = as<IntegerMatrix>(dat["capt_bin_unique"]);
  IntegerVector capt_bin_freqs = as<IntegerVector>(dat["capt_bin_freqs"]);
  int first_calls = as<int>(dat["first_calls"]);
  double cutoff = as<double>(dat["cutoff"]);
  NumericVector cutoff_vec(1);
  cutoff_vec(0) = cutoff;
  double lower_cutoff = as<double>(dat["lower_cutoff"]);
  NumericVector lower_cutoff_vec(1);
  lower_cutoff_vec(0) = lower_cutoff;
  double first_calls_trunc = as<double>(dat["first_calls_trunc"]);
  NumericMatrix capt_ss = as<NumericMatrix>(dat["capt_ss"]);
  NumericMatrix dists = as<NumericMatrix>(dat["dists"]);
  int trace = as<int>(dat["trace"]);
  double A = as<double>(dat["A"]);
  double dbl_min = as<double>(dat["DBL_MIN"]);
  double f = 0.0;
  int i;
  int j;
  double mu_ss;
  double det_prob;
  double sub_det_prob;
  double undet_prob;
  double undet_lower_prob;
  double sum_det_probs = 0;
  double sum_sub_det_probs = 0;
  double den;
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
      mu_ss = b0_ss - b1_ss*dists(i, j);
      expected_ss(j, i) = mu_ss;
      capt_prob = 1 - pnorm(cutoff_vec, mu_ss, sigma_ss);
      log_capt_probs(j, i) = log(capt_prob + dbl_min);
      log_evade_probs(j, i) = log(1 - capt_prob + dbl_min);
      undet_prob *= 1 - capt_prob;
      if (first_calls){
        undet_lower_prob *= pnorm(lower_cutoff_vec, mu_ss, sigma_ss);
      }
    }
    det_prob = 1 - undet_prob;
    sum_det_probs += det_prob;
    if (first_calls){
      // Saving detection probabilities for first calls, and
      // detection probabilities for any subsequent calls.
      mask_det_probs(i) = 0;
      mask_all_det_probs(i) = 0;
      if (1 - undet_lower_prob > first_calls_trunc){
	// Will probably find some numerical instability here.
	sum_sub_det_probs += (det_prob*undet_lower_prob)/(1 - undet_lower_prob);
	sub_det_prob = (det_prob*undet_lower_prob)/(1 - undet_lower_prob);
	mask_det_probs(i) += det_prob;
	mask_all_det_probs(i) += det_prob + sub_det_prob;
      }
    }
  }
  return D;
}

