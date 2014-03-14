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
//' dists <- distances(example.traps, example.mask)
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
//' bears <- bearings(example.traps, example.mask) 
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
NumericMatrix make_toa_ssq(const NumericMatrix& capt, const NumericMatrix& dists){
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
	  delts(index) = capt(i, k) - dists(k, j)/330;
	  index++;
	}
	out(i, j) = sum(pow(delts - mean(delts), 2));
      }
    }
  }  
  return out;
}
