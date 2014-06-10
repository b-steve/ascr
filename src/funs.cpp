#include <Rcpp.h>
using namespace Rcpp;

//' Calculating distances between two sets of points
//'
//' Calculates pairwise distances between two sets of points.
//'
//' @param a matrix containing a set of coordinates.
//' @param b matrix containing another set of coordinates.
//'
//' @return A matrix with pairwise distances between the two sets of points.
//'
//' @examples
//' dists <- distances(example$traps, example$mask)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix distances(const NumericMatrix& a, const NumericMatrix& b){
  int n_a = a.nrow();
  int n_b = b.nrow();
  NumericMatrix out(n_a, n_b);
  for (int j = 0; j < n_b; j++){
    for (int i = 0; i < n_a; i++){
      out(i, j) = pow(pow(a(i, 0) - b(j, 0), 2) + pow(a(i, 1) - b(j, 1), 2), 0.5);
    }
  }
  return out;
}

//' Calculating bearings between two sets of points
//'
//' Calculates the bearing from one set of points to another.
//'
//' @param a matrix containing a set of coordinates.
//' @param b matrix containing another set of coordinates.
//'
//' @return A matrix with bearings of the points in matrix b from the points in matrix a.
//'
//' @examples
//' bears <- bearings(example$traps, example$mask) 
//'
//' @export
// [[Rcpp::export]]
NumericMatrix bearings(const NumericMatrix& a, const NumericMatrix& b){
  int n_a = a.nrow();
  int n_b = b.nrow();
  NumericMatrix out(n_a, n_b);
  double x_diff, y_diff, bearing;
  for (int j = 0; j < n_b; j++){
    for (int i = 0; i < n_a; i++){
      x_diff = b(j, 0) - a(i, 0);
      y_diff = b(j, 1) - a(i, 1);
      bearing = atan(x_diff/y_diff);
      if (y_diff < 0){
	bearing += M_PI;
      } else if (x_diff < 0){
	bearing += 2*M_PI;
      }
      out(i, j) = bearing;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix make_toa_ssq(const NumericMatrix& capt, const NumericMatrix& dists, const double& sound_speed){
  int n = capt.nrow();
  int n_traps = capt.ncol();
  int n_mask = dists.ncol();
  NumericMatrix out(n, n_mask);
  int n_dets;
  int index;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n_mask; j++){
      n_dets = 0;
      for (int k = 0; k < n_traps; k++){
	if (capt(i, k) > 0) n_dets++;
      }
      NumericVector delts(n_dets);
      index = 0;
      for (int k = 0; k < n_traps; k++){
	if (capt(i, k) > 0){
	  delts(index) = capt(i, k) - dists(k, j)/sound_speed;
	  index++;
	}
	out(i, j) = sum(pow(delts - mean(delts), 2));
      }
    }
  }  
  return out;
}

// [[Rcpp::export]]
List find_local(const IntegerMatrix& capt, const NumericMatrix& dists, const double& buffer){
  int n = capt.nrow();
  int n_traps = capt.ncol();
  int n_mask = dists.ncol();
  //IntegerVector capt_hist(n_traps);
  IntegerVector capt_hist(n_traps);
  IntegerVector local_points(n_mask);
  int is_local;
  int n_local;
  int i, j, k;
  List out(n);
  for (i = 0; i < n; i++){
    // Getting current capture history.
    for (j = 0; j < n_traps; j++){
      capt_hist(j) = capt(i, j);
    }
    for (j = 0; j < n_mask; j++){
      // Resetting local_points to 0.
      local_points(j) = 0;
      is_local = 1;
      k = 0;
      while ((k < n_traps) & (is_local == 1)){
	// If a trap at which detection was made is not local to mask
	// point, set is_local to false.
	if (capt_hist(k) == 1){
	  if (dists(k, j) > buffer){
	    is_local = 0;
	  }
	}
	k++;
      }
      // Set locality of jth point.
      local_points(j) = is_local;
    }
    // Number of points local to ith detection.
    n_local = std::accumulate(local_points.begin(), local_points.end(), 0.0);
    IntegerVector which_local(n_local);
    k = 0;
    for (j = 0; j < n_mask; j++){
      if (local_points(j) == 1){
	which_local(k) = j + 1;
	k++;
      }
    }
    out[i] = which_local;
  }
  return out;
}
