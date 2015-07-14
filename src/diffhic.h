#ifndef HICCUP_H
#define HICCUP_H

#include <deque>
#include <queue>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "R.h"
#include "Rinternals.h"

template <class T>
struct sort_row_index {
	sort_row_index(const T* p) : ptr(p) {}
	bool operator() (const int& l, const int& r) const { 
		if (ptr[l]==ptr[r]) { return (l < r); }
		else { return (ptr[l] < ptr[r]); }
	}
private:
	const T* ptr;
};

extern "C" {

SEXP check_input(SEXP, SEXP);

SEXP cap_input(SEXP, SEXP, SEXP);

SEXP cluster_2d (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

SEXP split_clusters (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_bounding_box (SEXP, SEXP, SEXP);

SEXP quadrant_bg (SEXP, SEXP, 
	SEXP, SEXP, SEXP, 
	SEXP, SEXP, 
	SEXP, SEXP, SEXP);

SEXP count_background(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, 
		SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP count_connect(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP count_patch(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP iterative_correction(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP report_hic_pairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

SEXP report_hic_binned_pairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

SEXP test_parse_cigar(SEXP, SEXP);

SEXP test_fragment_assign(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP pair_stats (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_missing_dist(SEXP, SEXP, SEXP, SEXP);

SEXP directionality(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif
