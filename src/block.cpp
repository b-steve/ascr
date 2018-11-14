#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// Hard coding sound speed limits.
const double lower_limit = 150;
const double mid_limit = 310;
const double sound_speed = 330;

// Some prototypes.
int min_unallocated(const LogicalVector& allocated);
void add_to_block(const int& i, LogicalVector& in_block, const LogicalMatrix& mat);
void reset_in_block(LogicalVector& in_block);
NumericVector which_max_reqss(const NumericMatrix& reqss, const LogicalMatrix& allocated, IntegerMatrix& skip);
int min_skip_matrix(const IntegerMatrix& mat, const LogicalMatrix& allocated);
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

void set_false(const int& a, const int& b, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, bool& abort, const int& orig_a, const int& orig_b, const int& skip_val, IntegerMatrix& score_true, IntegerMatrix& score_false){
  int n = out.nrow();;
  // Abort if setting to false when > sound_speed.
  if (reqss(a, b) > sound_speed){
    abort = true;
    if (skip_val == 2){
      score_false(orig_a, orig_b) += 1;
      score_false(orig_b, orig_a) += 1;
    }
  } else {
    if (skip_val == 2){
      score_true(orig_a, orig_b) += 1;
      score_true(orig_b, orig_a) += 1;
    }
  }
  // Set to false.
  set_value(a, b, false, out, allocated);
  int i, k, fixed, unfixed;
  // Need to set all matches with each to false.
  for (k = 0; k < 2; k++){
    if (k == 0){
      fixed = a;
      unfixed = b;
    } else {
      fixed = b;
      unfixed = a;
    }
    for (i = 0; i < n; i++){
      if (allocated(i, unfixed) & !allocated(fixed, i)){
	// Checking if match.
	if (out(i, unfixed)){
	  // If so, set all comparisons between the matches to false.
	  set_false(fixed, i, out, allocated, reqss, abort, orig_a, orig_b, skip_val, score_true, score_false);
	}
      }
    }
  }
}

void set_true(const int& a, const int& b, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, bool& abort, const int& orig_a, const int& orig_b, const int& skip_val, IntegerMatrix& score_true, IntegerMatrix& score_false){
  int n = out.nrow();
  double limit;
  // Limit is sound_speed at first stage, mid_limit at second stage.
  if (skip_val == 0){
    limit = sound_speed;
  } else if (skip_val >= 1){
    limit = mid_limit;
  }
  // Abort if setting to true when < sound_speed (first
  // stage), or if setting to true when < mid_limit (second
  // stage).
  if (reqss(a, b) < limit){
    abort = true;
    if (skip_val == 2){
      score_false(orig_a, orig_b) += 1;
      score_false(orig_b, orig_a) += 1;
    }
  } else {
    if (skip_val == 2){
      score_true(orig_a, orig_b) += 1;
      score_true(orig_b, orig_a) += 1;
    }
  }
  // Set to true.
  set_value(a, b, true, out, allocated);
  int i, k, fixed, unfixed;
  // Need all other comparisons to match.
  for (k = 0; k < 2; k++){
    if (k == 0){
      fixed = a;
      unfixed = b;
    } else {
      fixed = b;
      unfixed = a;
    }
    for (i = 0; i < n; i++){
      // Match grouping elements between two cues.
      if (allocated(i, unfixed) & !allocated(fixed, i)){
	if (out(i, unfixed)){
	  set_true(fixed, i, out, allocated, reqss, abort, orig_a, orig_b, skip_val, score_true, score_false);
	} else {
	  set_false(fixed, i, out, allocated, reqss, abort, orig_a, orig_b, skip_val, score_true, score_false);
	}
      }
    }
  }
}

// Considers setting candidate element to true. Basically creates a copy of things and sees if things play out.
void consider(const int& a, const int& b, LogicalMatrix& out, LogicalMatrix& allocated, const NumericMatrix& reqss, IntegerMatrix& skip, IntegerMatrix& score_true, IntegerMatrix& score_false){
  int n = out.nrow();
  LogicalMatrix out_copy(n, n);
  LogicalMatrix allocated_copy(n, n);
  copy_matrix(out, out_copy);
  copy_matrix(allocated, allocated_copy);
  bool abort = false;
  // If up to fourth skip make a decision based on true/false scores.
  if ((skip(a, b) == 4) & (score_true(a, b) <= score_false(a, b))){
    set_false(a, b, out_copy, allocated_copy, reqss, abort, a, b, skip(a, b), score_true, score_false);
    abort = false;
  } else {
    // Otherwise do as usual.
    set_true(a, b, out_copy, allocated_copy, reqss, abort, a, b, skip(a, b), score_true, score_false);
    // If up to third skip, overwrite abort if winning on score.
    if ((skip(a, b) == 3) & (score_true(a, b) > score_false(a, b))){
      abort = false;
    }
  }
  if (!abort){
    copy_matrix(out_copy, out);
    copy_matrix(allocated_copy, allocated);
    reset_matrix(skip);
    reset_matrix(score_true);
    reset_matrix(score_false);
  } else {
    skip(a, b) = skip(a, b) + 1;
    skip(b, a) = skip(b, a) + 1;
  }
}

// Algorithm to resolve incomplete blocks.
// [[Rcpp::export]]
LogicalMatrix blockify(const LogicalMatrix& block, const NumericMatrix& reqss){
  int n = block.nrow();
  LogicalMatrix out(n, n);
  LogicalMatrix allocated(n, n);
  IntegerMatrix skip(n, n);
  IntegerMatrix score_true(n, n);
  IntegerMatrix score_false(n, n);
  NumericVector candidate(2);
  int i, j;
  // Allocating diagonal to TRUE and reqss < lower_limit as FALSE.
  for (i = 0; i < n; i++){
    for (j = i; j < n; j++){
      if (i == j){
	set_value(i, j, true, out, allocated);
      } else {
	if (reqss(i, j) < lower_limit){
	  set_value(i, j, false, out, allocated);
	} else {
	  out(i, j) = NA_LOGICAL;
	  out(j, i) = NA_LOGICAL;
	}
      }
    }
  }
  bool cont = true;
  while(cont){
    candidate = which_max_reqss(reqss, allocated, skip);
    int a, b;
    a = candidate(0);
    b = candidate(1);
    consider(a, b, out, allocated, reqss, skip, score_true, score_false);
    cont = is_false(all(allocated));
  }
  return out;
}

// Finds the earliest unallocated cue.
int min_unallocated(const LogicalVector& allocated){
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
    if ((!in_block(j)) & mat(i, j)){
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
  int min_skip = min_skip_matrix(skip, allocated);
  double current_max = 0;
  // If not up to fourth skip, choose largest reqss.
  for (i = 0; i < n - 1; i++){
    for (j = i + 1; j < n; j++){
      if ((reqss(i, j) >= current_max) & (skip(i, j) == min_skip) & (!allocated(i, j))){
	out(0) = i;
	out(1) = j;
	current_max = reqss(i, j);
      }
    }
  }
  // Otherwise choose smallest largest reqss.
  if (min_skip == 4){
    double current_min = current_max;
    for (i = 0; i < n - 1; i++){
      for (j = i + 1; j < n; j++){
	if ((reqss(i, j) <= current_min) & (skip(i, j) == min_skip) & !allocated(i, j)){
	  out(0) = i;
	  out(1) = j;
	  current_min = reqss(i, j);
	}
      }
    }
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
int min_skip_matrix(const IntegerMatrix& skip, const LogicalMatrix& allocated){
  int n = skip.nrow();
  int i, j;
  int out = 0;
  bool initial_found = false;
  for (i = 0; i < n - 1; i++){
    for (j = i + 1; j < n; j++){
      if (initial_found){
	if ((skip(i, j) < out) & !allocated(i, j)){
	  out = skip(i, j);
	}
      } else {
	if (!allocated(i, j)){
	  out = skip(i, j);
	  initial_found = true;
	}
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
