#include "diffhic.h"

struct orderer {
	orderer (const double* x) : cptr(x) {}
	bool operator() (const int& l, const int& r) const {  return cptr[l] < cptr[r]; }
	const double* cptr;
};

SEXP iterative_correction(SEXP avecount, SEXP anchor, SEXP target, SEXP local, 
		SEXP nfrags, SEXP iter, SEXP exlocal, SEXP lowdiscard, SEXP winsorhigh) try {

	// Checking vector type and length.
	if (!isNumeric(avecount)) { throw std::runtime_error("average counts must be supplied as a double-precision vector"); }
	if (!isInteger(anchor) || !isInteger(target)) { throw std::runtime_error("anchor/target indices must be supplied as an integer vector"); }
	if (!isLogical(local)) { throw std::runtime_error("local specifier should be a logical vector"); }
	const int npairs=LENGTH(avecount);
	if (npairs!=LENGTH(anchor) || npairs!=LENGTH(target) || npairs!=LENGTH(local)) { throw std::runtime_error("lengths of vectors are not equal"); }
	const double* acptr=REAL(avecount);
	const int* aptr=INTEGER(anchor),
		  * tptr=INTEGER(target),
		  * lptr=LOGICAL(local);

	// Checking scalar type.
	if (!isInteger(nfrags) || LENGTH(nfrags)!=1) { throw std::runtime_error("number of fragments should be an integer scalar"); }
	if (!isInteger(iter) || LENGTH(iter)!=1) { throw std::runtime_error("number of iterations should be an integer scalar"); }
	if (!isInteger(exlocal) || LENGTH(exlocal)!=1) { throw std::runtime_error("exclusion specifier should be an integer scalar"); }
	const int numfrags=asInteger(nfrags),
		  iterations=asInteger(iter),
		  excluded=asInteger(exlocal);
	if (!isNumeric(lowdiscard) || LENGTH(lowdiscard)!=1) { throw std::runtime_error("proportion to discard should be a double-precision scalar"); }
	const double discarded=asReal(lowdiscard);
	if (!isNumeric(winsorhigh) || LENGTH(winsorhigh)!=1) { throw std::runtime_error("proportion to winsorize should be a double-precision scalar"); }
	const double winsorized=asReal(winsorhigh);
	
	SEXP output=PROTECT(allocVector(VECSXP, 3));
try {
	// Setting up a SEXP to hold the working probabilities (ultimately the truth).
	SET_VECTOR_ELT(output, 0, allocVector(REALSXP, npairs));
	double* wptr=REAL(VECTOR_ELT(output, 0));
	std::copy(acptr, acptr+npairs, wptr);

	// Excluding highly local interactions.
	int num_values=npairs;
	bool drop_intras=(excluded==NA_INTEGER);
	if (drop_intras || excluded >= 0) {
		for (int j=0; j<npairs; ++j) {
			if (lptr[j] && (drop_intras || aptr[j]-tptr[j] <= excluded)) { 
				wptr[j]=R_NaReal; 
				--num_values;
			}
		}
	}

	/**************************************************
 	 * Getting rid of unstable or extreme elements.
 	 **************************************************/

	// Winsorizing for the most abundant interactions.
	const int todrop=int(double(num_values)*winsorized);
	if (todrop>0) {
		int* ordered=new int[npairs];
		try {
			for (int i=0; i<npairs; ++i) { ordered[i]=i; }
			orderer current(acptr);
			std::sort(ordered, ordered+npairs, current);
	
			// Identifying the winsorizing value.		
			int curcount=0;
 		   	double winsor_val=R_NaReal;
			for (int i=npairs-1; i>=0; --i) { 
				const int& curo=ordered[i];
				if (ISNA(wptr[curo])) { continue; }
				++curcount;
				if (curcount > todrop) { 
					winsor_val=wptr[curo];
					break; 
				}
			}

			// Applying said value to all high-abundance interactions.
			if (ISNA(winsor_val)) { throw std::runtime_error("specified winsorizing proportion censors all data"); }
			curcount=0;
			for (int i=npairs-1; i>=0; --i) { 
				const int& curo=ordered[i];
				if (ISNA(wptr[curo])) { continue; }
				wptr[curo]=winsor_val;
				++curcount;
				if (curcount > todrop) { break; }
			}		
		} catch (std::exception& e){
			delete [] ordered;
			throw;
		}
		delete[] ordered;
	}

	// Bias and coverage vectors.
	SET_VECTOR_ELT(output, 1, allocVector(REALSXP, numfrags));
	double* bias=REAL(VECTOR_ELT(output, 1)),
		* biaptr=bias-1;
	for (int i=0; i<numfrags; ++i) { bias[i]=R_NaReal; }
	for (int pr=0; pr<npairs; ++pr) { 
		if (!ISNA(wptr[pr])) { biaptr[aptr[pr]]=biaptr[tptr[pr]]=1; }
	}
	double* coverage=(double*)R_alloc(numfrags, sizeof(double));
	for (int i=0; i<numfrags; ++i) { coverage[i]=0; }
	double* covptr=coverage-1; // To deal with 1-based indices.

	// Removing low-abundance fragments.
	if (discarded > 0) {
		for (int pr=0; pr<npairs; ++pr) {  // Computing the coverage.
			if (ISNA(wptr[pr])) { continue; }
			const int& cura=aptr[pr];
			const int& curt=tptr[pr];
			covptr[cura]+=wptr[pr]; 
			if (cura!=curt) { covptr[curt]+=wptr[pr]; }
        }

		int* ordering=new int [numfrags];
		try {
			for (int i=0; i<numfrags; ++i) { ordering[i]=i; }
			orderer current(coverage);
			std::sort(ordering, ordering+numfrags, current);
		
			// Ignoring empty rows.
			int counter=0;
			while (counter < numfrags && ISNA(bias[ordering[counter]])){ ++counter; }
				
			// Filtering out low-abundance rows.
			const int leftover=int(discarded*double(numfrags-counter))+counter;
			while (counter < leftover) { 
				bias[ordering[counter]]=R_NaReal;
				++counter; 
			}

			// Setting everything else back to zero.
			while (counter < numfrags) { 
				coverage[ordering[counter]]=0;
				++counter;
			}
		} catch (std::exception& e) { 
			delete [] ordering;
			throw;
		}
		delete [] ordering;
	
		// Propogating the filters through the interactions to remove anything with an anchor in the removed fragments.
		for (int pr=0; pr<npairs; ++pr) {
			if (ISNA(wptr[pr])) { continue; }
			if (ISNA(biaptr[aptr[pr]]) || ISNA(biaptr[tptr[pr]])) { wptr[pr]=R_NaReal; }  
		}
		for (int i=0; i<numfrags; ++i) { bias[i]=R_NaReal; }
		for (int pr=0; pr<npairs; ++pr) { 
			if (!ISNA(wptr[pr])) { biaptr[aptr[pr]]=biaptr[tptr[pr]]=1; }
		}
	}
	
	// Something to hold a diagnostic (i.e., the maximum step).
	SET_VECTOR_ELT(output, 2, allocVector(REALSXP, iterations));
	double* optr=REAL(VECTOR_ELT(output, 2));

    /********************************************************
 	 * Now, actually performing the iterative correction. ***
 	 ********************************************************/

	for (int it=0; it<iterations; ++it) {
		for (int pr=0; pr<npairs; ++pr) {  // Computing the coverage (ignoring locals, if necessary).
			if (ISNA(wptr[pr])) { continue; }
			const int& cura=aptr[pr];
			const int& curt=tptr[pr];
			covptr[cura]+=wptr[pr]; 
			if (cura!=curt) { covptr[curt]+=wptr[pr]; }
        }
		
		/* Computing the 'aditional bias' vector, to take the mean across all fragments.
		 * The exact value of the mean doesn't really matter, it just serves to slow down 
		 * the algorithm to avoid overshooting the mark and lack of convergence. Everything
		 * ends up being the same so long as the column sums approach 1.
		 */
//		double mean_of_all=0;
//		for (int i=0; i<numfrags; ++i) { if (!ISNA(coverage[i])) { mean_of_all+=coverage[i]; } }
//		mean_of_all/=numfrags;
//		for (int i=0; i<numfrags; ++i) { if (!ISNA(coverage[i])) { coverage[i]/=mean_of_all; } }

		/* Okay, that didn't really work. The strategy that is described in the paper fails to
 		 * give me contact probabilities that sum to 1. So instead, I'm using the coverage
 		 * itself to divide the working probabilities. Square rooting is necessary to avoid
 		 * instability during iteration.
 		 */
		for (int i=0; i<numfrags; ++i) { if (!ISNA(bias[i])) { coverage[i]=std::sqrt(coverage[i]); } }

		// Dividing the working matrix with the (geometric mean of the) additional biases.	
		for (int pr=0; pr<npairs; ++pr){ 
			if (!ISNA(wptr[pr])) { wptr[pr]/=covptr[aptr[pr]]*covptr[tptr[pr]]; }
        }
		
		// Multiplying the biases by additional biases, and cleaning out coverage. We store
		// the maximum step to see how far off convergence it is.
		double& maxed=(optr[it]=0);
		for (int i=0; i<numfrags; ++i) {			
			if (!ISNA(bias[i])) {
				double& cur_cov=coverage[i];
				if (cur_cov > maxed) { maxed=cur_cov; }
				else if (1/cur_cov > maxed) { maxed=1/cur_cov; }
				bias[i]*=cur_cov;
				cur_cov=0;
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

