#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' Calculating distances between two sets of points
//'
//' Calculates pairwise distances between two sets of points.
//'
//' @param a matrix containing a set of coordinates.
//' @param b matrix containing another set of coordinates.
//' @return A matrix with pairwise distances between the two sets of points.
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
//' @return A matrix with bearings of the points in matrix b from the points in matrix a.
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
