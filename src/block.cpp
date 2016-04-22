#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// Some prototypes.
int min_unallocated(const LogicalVector& allocated);
void add_to_block(const int& i, LogicalVector& in_block, const LogicalMatrix& mat);
void reset_in_block(LogicalVector& in_block);
NumericVector which_max_reqss(const NumericMatrix& reqss, const LogicalMatrix& allocated);
void consider_set_true(const NumericVector& candidate, const LogicalMatrix& out, const LogicalMatrix& allocated, const NumericMatrix& reqss);
void set_true(const NumericVector& candidate, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, bool& abort);
void set_value(const NumericVector& i, const NumericVector& j, const bool& value, LogicalMatrix& out, LogicalMatrix& allocated);
  
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

// Algorithm to resolve incomplete blocks.
// [[Rcpp::export]]
LogicalMatrix blockify(const LogicalMatrix& block, const NumericMatrix& reqss){
  int n = block.nrow();
  LogicalMatrix out(n, n);
  LogicalMatrix allocated(n, n);
  NumericVector candidate(2);
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
	out(i, j) = true;
	allocated(i, j) = true;
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
  candidate = which_max_reqss(reqss, allocated);
  consider_set_true(candidate, out, allocated, reqss);
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
NumericVector which_max_reqss(const NumericMatrix& reqss, const LogicalMatrix& allocated){
  int n = reqss.nrow();
  NumericVector out(2);
  double current_max = 0;
  int i, j;
  for (i = 0; i < n - 1; i++){
    for (j = i + 1; j < n; j++){
      if (reqss(i, j) > current_max & !allocated(i, j)){
	out(0) = i;
	out(1) = j;
	current_max = reqss(i, j);
      }
    }
  }
  return out;
}

// Considers setting candidate element to true. Basically creates a copy of things and sees if things play out.
void consider_set_true(const NumericVector& candidate, const LogicalMatrix& out, const LogicalMatrix& allocated, const NumericMatrix& reqss){
  LogicalMatrix out_copy = out;
  LogicalMatrix allocated_copy = allocated;
  bool abort = false;
  set_true(candidate, out_copy, allocated_copy, reqss, abort);
}

void set_true(const NumericVector& candidate, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, bool& abort){
  int a = candidate(0);
  int b = candidate(1);
  out(a, b) = true;
  out(b, a) = true;
  allocated(a, b) = true;
  allocated(b, a) = true;
}

void set_value(const int& i, const int& j, const bool& value, LogicalMatrix& out, LogicalMatrix& allocated){
  out(i, j) = value;
  allocated(i, j) = true;
  if (i != j){
    out(j, i) = value;
    allocated(j, i) = true;
  }
}
