#include "diffhic.h"

/* This just provides a function that quickly and efficiently
 * checks the incoming pair counts for (a) anchor >= target,
 * (b) ordering of anchor and targets.
 */

SEXP check_input (SEXP anchor, SEXP target) try { 
	if (!isInteger(anchor)) { throw std::runtime_error("anchor should be an integer vector"); }
	if (!isInteger(target)) { throw std::runtime_error("target should be an integer vector"); }
	const int nlen=LENGTH(anchor);
	if (LENGTH(target)!=nlen) { throw std::runtime_error("vectors should be of the same length"); }

	const int * aptr=INTEGER(anchor),
		  * tptr=INTEGER(target);

	if (nlen) {
		if (aptr[0] < tptr[0]) { throw std::runtime_error("anchor should be greater than or equal to target"); }
		for (int i=1; i<nlen; ++i) {
			if (aptr[i] < tptr[i]) { throw std::runtime_error("anchor should be greater than or equal to target"); }
			if (aptr[i] < aptr[i-1] || (aptr[i]==aptr[i-1] && tptr[i] < tptr[i-1])) {
				throw std::runtime_error("pairs should be sorted by anchor and target"); 
			}
		}
	}
	return ScalarLogical(1);
} catch (std::exception& e){
	return mkString(e.what());
}

/* This provides a function that attempts to cap the 
 * number of read pairs generated from each restriction fragment.
 * It assumes that all anchor and target indices are already
 * sorted, based on passing check_input.
 */

SEXP cap_input (SEXP anchor, SEXP target, SEXP cap) try { 
	if (!isInteger(anchor)) { throw std::runtime_error("anchor should be an integer vector"); }
	if (!isInteger(target)) { throw std::runtime_error("target should be an integer vector"); }
	const int nlen=LENGTH(anchor);
	if (LENGTH(target)!=nlen) { throw std::runtime_error("vectors should be of the same length"); }

	const int * aptr=INTEGER(anchor),
		  * tptr=INTEGER(target);
	if (!isInteger(cap) || LENGTH(cap)!=1) { throw std::runtime_error("cap should be an integer scalar"); }
	const int recap=asInteger(cap);

	SEXP output=PROTECT(allocVector(LGLSXP, nlen));
try {
	if (nlen) { 
		int * optr=LOGICAL(output);
		optr[0]=1;
		int counter=1;
		for (int i=1; i<nlen; ++i) {
			if (aptr[i]==aptr[i-1] && tptr[i]==tptr[i-1]) {
				++counter;
				optr[i]=(counter <= recap);
			} else {
				counter=1;
				optr[i]=1;
			}
		}
	}
} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e){
	return mkString(e.what());
}
