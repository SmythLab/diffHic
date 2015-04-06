#include "diffhic.h"

extern "C" {

SEXP get_missing_dist (SEXP chrends, SEXP existing_anchor, SEXP existing_target, SEXP middies) try {
	if (!isInteger(chrends)) { throw std::runtime_error("chromosome end indices must be integer"); }
	const int nchrs = LENGTH(chrends);
	if (!isInteger(existing_anchor)) { throw std::runtime_error("anchor indices must be integer"); }
	if (!isInteger(existing_target)) { throw std::runtime_error("target indices must be integer"); }
	const int npts = LENGTH(existing_anchor);
	if (npts!=LENGTH(existing_target)) { throw std::runtime_error("anchor and target index vectors must be of the same length"); }
	if (!isNumeric(middies)) { throw std::runtime_error("midpoint vector must be double-precision"); }

	const int* cptr=INTEGER(chrends);
	const int* aptr=INTEGER(existing_anchor);
	const int* tptr=INTEGER(existing_target);
	const double* mptr=REAL(middies);
	
	std::deque<double> stored;
	int chr_start=0, chr_index=0, pt_index=0;
	int anchor=0, target=0;
	bool ispresent;

	while (chr_index < nchrs) {
 	   	const int& chr_end=cptr[chr_index];
		
		for (anchor=chr_start; anchor<chr_end; ++anchor) {
			for (target=chr_start; target<=anchor; ++target) {
				ispresent=false;
 			   	while (pt_index < npts && aptr[pt_index]==anchor && tptr[pt_index]==target) { 
					ispresent=true;
					++pt_index; 
				}
				if (!ispresent) { stored.push_back(mptr[anchor]-mptr[target]); }
			}
		}

		chr_start=chr_end;
		++chr_index;
	}
	if (pt_index!=npts) { throw std::runtime_error("failed to parse all supplied points"); }

	SEXP output=PROTECT(allocVector(REALSXP, stored.size()));
try{
	double * optr=REAL(output);
	for (size_t o=0; o<stored.size(); ++o) { optr[o]=stored[o]; }
} catch (std::exception& e) { 
	UNPROTECT(1);
	throw;
} 
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

}
