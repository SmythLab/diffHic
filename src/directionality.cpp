#include "read_count.h"

SEXP directionality(SEXP all, SEXP bin, SEXP span, SEXP first_bin, SEXP last_bin, 
		SEXP max_it, SEXP tolerance, SEXP offsets, SEXP dispersion) try {

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

	// Setting up the NB average stuff.
	if (!isInteger(max_it) || LENGTH(max_it)!=1) { throw std::runtime_error("maximum number of iterations must be an integer scalar"); }
	const int maxit=asInteger(max_it);
	if (!isReal(tolerance) || LENGTH(tolerance)!=1) { throw std::runtime_error("tolerance must be a double-precision vector"); }
	const double tol=asReal(tolerance);
	if (!isReal(offsets) || LENGTH(offsets)!=nlibs) { throw std::runtime_error("offsets must be a double-precision vector of length equal to number of libraries"); }
	const double* offptr=REAL(offsets);
	if (!isReal(dispersion) || LENGTH(dispersion)!=1) { throw std::runtime_error("dispersion must be a double-precision vector"); }
	const double disp=asReal(dispersion);

	// Setting up the memory containers.
	const int nbins=lbin-fbin+1;
	int* curcounts=(int*)R_alloc(nlibs*nbins, sizeof(int)); 
	bool* ischanged=(bool*)R_alloc(nbins, sizeof(bool));
	for (int i=0; i<nbins; ++i) { ischanged[i]=false; }
	std::deque<int> waschanged, counts, anchors, targets;

	// Setting up the output vectors immediately.
	SEXP output=PROTECT(allocVector(VECSXP, 2));
try {
	SET_VECTOR_ELT(output, 0, allocVector(REALSXP, nbins));
	double* downptr=REAL(VECTOR_ELT(output, 0));
	SET_VECTOR_ELT(output, 1, allocVector(REALSXP, nbins));
	double* upptr=REAL(VECTOR_ELT(output, 1));
	for (int i=0; i<nbins; ++i) { downptr[i]=upptr[i]=0; }

	// Other assorted sundries.
	size_t vecdex;
	int rowdex, curanchor;
	double current_average;
   	size_t diff;

	while (!engine.empty()) {
		engine.fill(curcounts, ischanged, waschanged);
		curanchor=engine.get_anchor() - fbin;
		for (vecdex=0; vecdex<waschanged.size(); ++vecdex) {
			rowdex=waschanged[vecdex];
	
			// Filling up the directionality indices.		
			diff=curanchor-rowdex;
			if (diff && diff <= sp) { 
				current_average=nb_average(nlibs, maxit, tol, offptr, curcounts+rowdex*nlibs, disp);
				downptr[curanchor]+=current_average;
				upptr[rowdex]+=current_average;
			} 

			// Resetting the ischanged vector for the next stretch of bin anchors.
			ischanged[waschanged[vecdex]]=false;
		}
		waschanged.clear();
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

