#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// Some prototypes.
int min_unallocated(const LogicalVector& allocated);
void add_to_block(const int& i, LogicalVector& in_block, const LogicalMatrix& mat);
void reset_in_block(LogicalVector& in_block);
NumericVector which_max_reqss(const NumericMatrix& reqss, const LogicalMatrix& allocated, IntegerMatrix& skip);
int min_matrix(const IntegerMatrix& mat);
void copy_matrix(const NumericMatrix& from, NumericMatrix& to);
void copy_matrix(const LogicalMatrix& from, LogicalMatrix& to);
void reset_matrix(IntegerMatrix& mat);

// Functions for determining cue sources.

// Finding sets of incomplete blocks.
// [[Rcpp::export]]
IntegerVector find_incomplete_blocks(const LogicalMatrix& mat){
  int n = mat.nrow();
  IntegerVector out(n);
  LogicalVector allocated(n);
  LogicalVector in_block(n);
  int i, j;
  int block = 1;
  bool cont = true;
  while (cont){
    // Resetting in_block.
    reset_in_block(in_block);
    // Taking first unallocated cue as a starting point.
    i = min_unallocated(allocated);
    // Adding all block members; alters in_block.
    add_to_block(i, in_block, mat);
    for (j = 0; j < n; j++){
      if (in_block(j)){
	out(j) = block;
	allocated(j) = true;
      }
    }
    block++;
    cont = is_false(all(allocated));
  }
  return out;
}

void set_value(const int& i, const int& j, const bool& value, LogicalMatrix& out, LogicalMatrix& allocated){
  out(i, j) = value;
  allocated(i, j) = true;
  if (i != j){
    out(j, i) = value;
    allocated(j, i) = true;
  }
}

void set_true(const int& a, const int& b, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, bool& abort){
  int n = out.nrow();
  // Set candidate to true.
  set_value(a, b, true, out, allocated);
  int i, k, fixed, unfixed;
  for (k = 0; k < 2; k++){
    if (k == 0){
      fixed = a;
      unfixed = b;
    } else {
      fixed = b;
      unfixed = a;
    }
    for (i = 0; i < n; i++){
      // Figure out which other elements associated with each a and b must be true.
      if (allocated(i, unfixed) & !allocated(fixed, i)){
	if (out(i, unfixed)){
	  // Abort if setting to true when < 330.
	  if (reqss(fixed, i) < 330){
	    abort = true;
	  }
	  set_true(fixed, i, out, allocated, reqss, abort);
	} else {
	  // Abort if setting to false when > 330.
	  if (abort = reqss(fixed, i) > 330){
	    abort = true;
	  }
	  set_value(fixed, i, false, out, allocated);
	}
      }
    }
  }
}

// Considers setting candidate element to true. Basically creates a copy of things and sees if things play out.
void consider_set_true(const int& a, const int& b, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, IntegerMatrix& skip){
  int n = out.nrow();
  LogicalMatrix out_copy(n, n);
  LogicalMatrix allocated_copy(n, n);
  copy_matrix(out, out_copy);
  copy_matrix(allocated, allocated_copy);
  bool abort = false;
  set_true(a, b, out_copy, allocated_copy, reqss, abort);
  if (!abort){
    copy_matrix(out_copy, out);
    copy_matrix(allocated_copy, allocated);
    reset_matrix(skip);
  } else {
    skip(a, b) = skip(a, b) + 1;
  }
}


// Algorithm to resolve incomplete blocks.
// [[Rcpp::export]]
LogicalMatrix blockify(const LogicalMatrix& block, const NumericMatrix& reqss){
  int n = block.nrow();
  LogicalMatrix out(n, n);
  LogicalMatrix allocated(n, n);
  IntegerMatrix skip(n, n);
  NumericVector candidate(2);
  bool set_to;
  int i, j;
  // Hard-coded absolute minimum speed.
  double min_ss = 275;
  // Hard-coded lower acceptable speed.
  double acc_ss = 310;
  // Hard-coded ideal acceptable speed.
  double ideal_ss = 330;
  // Allocating diagonal to TRUE and reqss < min_ss as FALSE.
  for (i = 0; i < n; i++){
    for (j = i; j < n; j++){
      if (i == j){
	set_value(i, j, true, out, allocated);
      } else {
	if (reqss(i, j) < min_ss){
	  set_value(i, j, false, out, allocated);
	} else {
	  out(i, j) = NA_LOGICAL;
	  out(j, i) = NA_LOGICAL;
	}
      }
    }
  }
  bool cont = true;
  cout << "Out:" << endl << out << endl;
  cout << "Allocated:" << endl << allocated << endl;
  cout << "Skip:" << endl << skip << endl;
  while(cont){
    candidate = which_max_reqss(reqss, allocated, skip);
    int a, b;
    a = candidate(0);
    b = candidate(1);
    consider_set_true(a, b, out, allocated, reqss, skip);
    cont = is_false(all(allocated));
    cout << "Candidate: " << candidate << endl;
    cout << "Out:" << endl << out << endl;
    cout << "Allocated:" << endl << allocated << endl;
    cout << "Skip:" << endl << skip << endl;
  }
  return out;
}

// Finds the earliest unallocated cue.
int min_unallocated(const LogicalVector& allocated){
  int n = allocated.size();
  int i = -1;
  bool cont = true;
  while (cont){
    i++;
    cont = allocated(i);
  }
  return i;
}

// Determines all cues in the same block as the ith cue.
void add_to_block(const int& i, LogicalVector& in_block, const LogicalMatrix& mat){
  int n = mat.nrow();
  in_block(i) = true;
  for (int j = 0; j < n; j++){
    if (!in_block(j) & mat(i, j)){
      add_to_block(j, in_block, mat);
    }
  }
}

// Resets in_block vector to all false.
void reset_in_block(LogicalVector& in_block){
  int n = in_block.size();
  for (int i = 0; i < n; i++){
    in_block(i) = false;
  }
}

// Generates distances between all detections.
// [[Rcpp::export]]
NumericMatrix detection_dists(const NumericMatrix& trap_dists, const NumericVector& traps){
  int n = traps.size();
  NumericMatrix out(n, n);
  for (int i = 0; i < n; i++){
    for (int j = i; j < n; j++){
      out(i, j) = trap_dists(traps(i) - 1, traps(j) - 1);
      out(j, i) = trap_dists(traps(i) - 1, traps(j) - 1);
    }
  }
  return out;
}

// Generates time differences between all detections.
// [[Rcpp::export]]
NumericMatrix detection_timediffs(const NumericVector& times, const NumericVector& traps){
  int n = traps.size();
  NumericMatrix out(n, n);
  for (int i = 0; i < n; i++){
    for (int j = i; j < n; j++){
      out(i, j) = abs(times(i) - times(j));
      out(j, i) = abs(times(i) - times(j));
    }
  }
  return out;
}

// Finds unallocated element of reqss with maximum speed.
NumericVector which_max_reqss(const NumericMatrix& reqss, const LogicalMatrix& allocated, IntegerMatrix& skip){
  int n = reqss.nrow();
  NumericVector out(2);
  int i, j;
  double current_max = 330;
  bool complete = false;
  int min_skip = min_matrix(skip);
  for (i = 0; i < n - 1; i++){
    for (j = i + 1; j < n; j++){
      if (reqss(i, j) > current_max & skip(i, j) == min_skip & !allocated(i, j)){
	out(0) = i;
	out(1) = j;
	current_max = reqss(i, j);
	complete = true;
      }
    }
  }
  if (!complete){
    cout << "Error: no allocations left above 330 ms." << endl;
    exit(1234);
  }
  return out;
}

void copy_matrix(const NumericMatrix& from, NumericMatrix& to){
  int n = from.nrow();
  int i, j;
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      to(i, j) = from(i, j);
    }
  }
}

void copy_matrix(const LogicalMatrix& from, LogicalMatrix& to){
  int n = from.nrow();
  int i, j;
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      to(i, j) = from(i, j);
    }
  }
}

// [[Rcpp::export]]
int min_matrix(const IntegerMatrix& mat){
  int n = mat.nrow();
  int i, j;
  int out = mat(0, 0);
  for (i = 0; i < n - 1; i++){
   for (j = i + 1; j < n; j++){
      if (mat(i, j) < out){
	out = mat(i, j);
      }
    }
  }
  return out;
}

void reset_matrix(IntegerMatrix& mat){
  int n = mat.nrow();
  int i, j;
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      mat(i, j) = 0;
    }
  }
}
