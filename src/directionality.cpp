#include "read_count.h"

SEXP directionality(SEXP all, SEXP bin, SEXP span, SEXP first_bin, SEXP last_bin) try {

	if (!isInteger(span) || LENGTH(span)!=1) { throw std::runtime_error("span to compute directionality must be an integer scalar"); }
	const size_t sp=asInteger(span);

	// Getting the indices of the first and last bin on the chromosome.
	if (!isInteger(first_bin) || LENGTH(first_bin)!=1) { throw std::runtime_error("index of first bin must be an integer scalar"); }
	const int fbin=asInteger(first_bin);
	if (!isInteger(last_bin) || LENGTH(last_bin)!=1) { throw std::runtime_error("index of last bin must be an integer scalar"); }
	const int lbin=asInteger(last_bin);

	// Setting up the binning engine.
	binner engine(all, bin, fbin, lbin);
	const int nlibs=engine.get_nlibs();
	const int nbins=engine.get_nbins();

    // Settin gup the memory containers.
	std::deque<int> counts, anchors, targets;

	// Setting up the output vectors immediately.
	SEXP output=PROTECT(allocVector(VECSXP, 2));
try {
	SET_VECTOR_ELT(output, 0, allocMatrix(INTSXP, nbins, nlibs));
	int** downptrs=(int**)R_alloc(nlibs, sizeof(int*));
	SET_VECTOR_ELT(output, 1, allocMatrix(INTSXP, nbins, nlibs));
	int** upptrs=(int**)R_alloc(nlibs, sizeof(int*));
    if (nlibs) {
        downptrs[0]=INTEGER(VECTOR_ELT(output, 0));
        upptrs[0]=INTEGER(VECTOR_ELT(output, 1));
        for (int i=1; i<nlibs; ++i) {
            downptrs[i]=downptrs[i-1]+nbins;
            upptrs[i]=upptrs[i-1]+nbins;
        }
	    for (int i=0; i<nlibs; ++i) {
            std::fill(downptrs[i], downptrs[i]+nbins, 0);
            std::fill(upptrs[i], upptrs[i]+nbins, 0);
        }
    }

	// Other assorted sundries.
	size_t vecdex;
	int rowdex, curanchor;
	double current_average;
   	size_t diff;
    int lib;
    std::deque<int>::const_iterator ccIt, wcIt;

	while (!engine.empty()) {
		engine.fill();
		curanchor=engine.get_anchor() - fbin;
        const std::deque<int>& waschanged=engine.get_changed();
        const std::deque<int>& curcounts=engine.get_counts();

        for (wcIt=waschanged.begin(); wcIt!=waschanged.end(); ++wcIt) {
       	    rowdex=(*wcIt);
			diff=curanchor-rowdex;

			// Filling up the directionality indices.		
			if (diff && diff <= sp) { 
                ccIt=curcounts.begin() + rowdex*nlibs;
                for (lib=0; lib<nlibs; ++lib, ++ccIt) {
                    const int& thiscount=(*ccIt);
                    downptrs[lib][curanchor]+=thiscount;
                    upptrs[lib][rowdex]+=thiscount;
                }
			} 
		}
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

