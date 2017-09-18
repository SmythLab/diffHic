#include "diffhic.h"
#include "utils.h"

SEXP get_missing_dist (SEXP chrends, SEXP existing_anchor1, SEXP existing_anchor2, SEXP middies) {
    BEGIN_RCPP

    Rcpp::IntegerVector ce(chrends), ea1(existing_anchor1), ea2(existing_anchor2);
    Rcpp::NumericVector mid(middies);
	const int nchrs = ce.size();
	const int npts = ea1.size();
	if (npts!=ea2.size()) { throw std::runtime_error("first and second anchor index vectors must be of the same length"); }

	std::deque<double> stored;
	int chr_start=0, chr_index=0, pt_index=0;

	while (chr_index < nchrs) {
 	   	const int& chr_end=ce[chr_index];
		
		for (int anchor=chr_start; anchor<chr_end; ++anchor) {
			for (int target=chr_start; target<=anchor; ++target) {
				bool ispresent=false;
 			   	while (pt_index < npts && ea1[pt_index]==anchor && ea2[pt_index]==target) { 
					ispresent=true;
					++pt_index; 
				}
				if (!ispresent) { stored.push_back(mid[anchor]-mid[target]); }
			}
		}

		chr_start=chr_end;
		++chr_index;
	}
	if (pt_index!=npts) { throw std::runtime_error("failed to parse all supplied points"); }

    return Rcpp::NumericVector(stored.begin(), stored.end());
    END_RCPP
}
