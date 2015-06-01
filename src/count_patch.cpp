#include "read_count.h"

SEXP count_patch(SEXP all, SEXP bin, SEXP filter, SEXP firstbin, SEXP lastbin) try {
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int f=asInteger(filter);

	// Getting the indices of the first and last bin on the target chromosome.
	if (!isInteger(firstbin) || LENGTH(firstbin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int fbin=asInteger(firstbin);
	if (!isInteger(lastbin) || LENGTH(lastbin)!=1) { throw std::runtime_error("index of last bin on target chromosome must be an integer scalar"); }
	const int lbin=asInteger(lastbin);

	// Setting up the binning engine.
	binner engine(all, bin, fbin, lbin);
	const int nlibs=engine.get_nlibs();

	// Setting up the memory containers.
	const int nbins=lbin-fbin+1;
	int* curcounts=(int*)R_alloc(nlibs*nbins, sizeof(int)); 
	bool* ischanged=(bool*)R_alloc(nbins, sizeof(bool));
	for (int i=0; i<nbins; ++i) { ischanged[i]=false; }
	std::deque<int> waschanged, counts, anchors, targets;

	// Sundries.
	size_t vecdex;
	int rowdex, countsum, curlib, curanchor;

	// Running through all libraries.
	while (!engine.empty()) {
		engine.fill(curcounts, ischanged, waschanged);
		curanchor=engine.get_anchor();

		// Adding it to the main list, if it's large enough.
		for (vecdex=0; vecdex<waschanged.size(); ++vecdex) {
			rowdex=waschanged[vecdex]*nlibs;
			countsum=0;
			for (curlib=0; curlib<nlibs; ++curlib) { countsum+=curcounts[rowdex+curlib]; }
			if (countsum >= f) { 
				anchors.push_back(curanchor);
				targets.push_back(waschanged[vecdex]+fbin);
				for (curlib=0; curlib<nlibs; ++curlib) { counts.push_back(curcounts[rowdex+curlib]); }
			}

			// Resetting the ischanged vector for the next stretch of bin anchors.
			ischanged[waschanged[vecdex]]=false;
		}
		waschanged.clear();
	}

	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		const int ncombos=anchors.size();
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ncombos));
		int* aoptr=INTEGER(VECTOR_ELT(output, 0));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ncombos));
		int* toptr=INTEGER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, allocMatrix(INTSXP, ncombos, nlibs));
		std::deque<int*> coptrs(nlibs);
		for (curlib=0; curlib<nlibs; ++curlib) {
			if (curlib==0) { coptrs[curlib]=INTEGER(VECTOR_ELT(output, 2)); }
			else { coptrs[curlib]=coptrs[curlib-1]+ncombos; }
		}	
		
		// Iterating across and filling both the matrix and the components.
		int cdex=-1;
		for (vecdex=0; vecdex<anchors.size(); ++vecdex) {
			aoptr[vecdex]=anchors[vecdex];
			toptr[vecdex]=targets[vecdex];
			for (curlib=0; curlib<nlibs; ++curlib) { coptrs[curlib][vecdex]=counts[++cdex]; }
		}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
