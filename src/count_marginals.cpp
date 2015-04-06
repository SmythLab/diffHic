#include "diffhic.h"

SEXP count_marginals (SEXP all, SEXP bins, SEXP total) try {
	if (!isInteger(bins)) { throw std::runtime_error("bin indices must be an integer vector"); }
	const int* bptr=INTEGER(bins)-1;
	const int nb=LENGTH(bins);
   	if (!isNewList(all)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }
	const int nlibs=LENGTH(all);
	if (!isInteger(total)||LENGTH(total)!=1) { throw std::runtime_error("total number of fragments should be an integer scalar"); }
	const int ntot=asInteger(total);

	// Setting up output for storage.
	SEXP output=PROTECT(allocMatrix(INTSXP, ntot, nlibs));
	try {
		int** optrs=(int**)R_alloc(nlibs, sizeof(int*));
		optrs[0]=INTEGER(output)-1;
		for (int i=1; i<nlibs; ++i) { optrs[i]=optrs[i-1]+ntot; }
		for (int i=0; i<nlibs; ++i) { for (int j=1; j<=ntot; ++j) { optrs[i][j]=0; } }

    	for (int i=0; i<nlibs; ++i) {
        	SEXP current=VECTOR_ELT(all, i);
        	if (!isNewList(current) || LENGTH(current)!=2) { 
				throw std::runtime_error("interactions must be supplied as a data.frame with anchor.id, target.id"); }
			
			const int* aaptr, *ttptr;
			int num=0;
        	for (int j=0; j<2; ++j) {
            	SEXP current_col=VECTOR_ELT(current, j);
            	if (!isInteger(current_col)) { throw std::runtime_error("interaction data must be in integer format"); }
            	int* ptr=INTEGER(current_col);
            	switch (j) {
                	case 0: 
						aaptr=ptr; 
						num=LENGTH(current_col);
						break;
                	case 1: ttptr=ptr; break;
                	default: break;
            	}
			}

			// Running through each row and adding it to the sum for the fragment.
			for (int j=0; j<num; ++j) {
				const int& cur_a=aaptr[j];
				const int& cur_t=ttptr[j];
				if (cur_a > nb || cur_t > nb) { throw std::runtime_error("anchor or target indices out of range for conversion to bin index"); }
				const int& aid=bptr[cur_a];
				optrs[i][aid]+=1;
				const int& tid=bptr[cur_t];
				if (tid==aid) { continue; }
				optrs[i][tid]+=1;
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
