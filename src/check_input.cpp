#include "diffhic.h"
#include "utils.h"

/* This just provides a function that quickly and efficiently
 * checks the incoming pair counts for (a) anchor >= target,
 * (b) ordering of anchor and targets.
 */

SEXP check_input (SEXP anchor1, SEXP anchor2) { 
    BEGIN_RCPP

    const Rcpp::IntegerVector a1(anchor1), a2(anchor2);
	const int nlen=a1.size();
	if (a2.size()!=nlen) { 
        throw std::runtime_error("vectors should be of the same length"); 
    }

	if (nlen) {
		if (a1[0] < a2[0]) { throw std::runtime_error("anchor should be greater than or equal to target"); }
		for (int i=1; i<nlen; ++i) {
			if (a1[i] < a2[i]) { throw std::runtime_error("anchor should be greater than or equal to target"); }
			if (a1[i] < a1[i-1] || (a1[i]==a1[i-1] && a2[i] < a2[i-1])) {
				throw std::runtime_error("pairs should be sorted by anchor and target"); 
			}
		}
	}
    
    return Rcpp::LogicalVector::create(true);
    END_RCPP
}

/* This provides a function that attempts to cap the 
 * number of read pairs generated from each restriction fragment.
 * It assumes that all anchor and target indices are already
 * sorted, based on passing check_input.
 */

SEXP cap_input (SEXP anchor1, SEXP anchor2, SEXP cap) { 
    BEGIN_RCPP
    const Rcpp::IntegerVector a1(anchor1), a2(anchor2);
    const int nlen=a1.size();
    int recap=check_integer_scalar(cap, "cap");

    Rcpp::LogicalVector output(nlen, true);
	if (nlen) { 
		int counter=1;
		for (int i=1; i<nlen; ++i) {
			if (a1[i]==a1[i-1] && a2[i]==a2[i-1]) {
				++counter;
                if (counter > recap) {
                    output[i]=false;
                }
			} else {
				counter=1;
			}
		}
	}

	return output;
    END_RCPP
}
