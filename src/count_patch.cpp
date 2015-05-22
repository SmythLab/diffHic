#include "read_count.h"

binner::binner(SEXP all, SEXP bin, int f, int l) : fbin(f), lbin(l) {
	if (!isInteger(bin)) { throw std::runtime_error("anchor bin indices must be integer vectors"); }
	bptr=INTEGER(bin)-1; // Assuming 1-based indices for anchors and targets.

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
	if (!isNewList(all)) { throw std::runtime_error("data on interacting read pairs must be contained within a list"); }
	nlibs=LENGTH(all);
	aptrs.resize(nlibs);
	tptrs.resize(nlibs);
	indices.resize(nlibs);
	nums.resize(nlibs);
	
	for (int i=0; i<nlibs; ++i) {
		SEXP current=VECTOR_ELT(all, i);
		if (!isNewList(current) || LENGTH(current)!=2) { 
			throw std::runtime_error("interactions must be supplied as a data.frame with anchor.id and target.id"); }

		for (int j=0; j<2; ++j) {
			SEXP current_col=VECTOR_ELT(current, j);
			if (!isInteger(current_col)) { throw std::runtime_error("interaction data must be in integer format"); }
			const int* ptr=INTEGER(current_col);
			switch (j) {
				case 0: 
					aptrs[i]=ptr; 
					nums[i]=LENGTH(current_col);
					break;
				case 1: 
					tptrs[i]=ptr; 
					if (LENGTH(current_col)!=nums[i]) { throw std::runtime_error("vectors should be the same length"); }		
					break;
				default: break;
			}
		}
		
		// Populating the priority queue.
		if (nums[i]) { next.push(coord(bptr[aptrs[i][0]], bptr[tptrs[i][0]], i)); }
	}
	return;
}

void binner::fill(int* curcounts, bool* ischanged, std::deque<int>& waschanged) {
	/* Running through all libraries. The idea is to use stretches of identical bin anchors, such that
	 * we only need to worry about different bin targets (i.e., the problem becomes 1-dimensional).
	 * This assumes that the bin transformation is monotonic, and that anchors are sorted. The function 
	 * will stop once all identical bin anchors have been processed. Counts for all bins 
	 * in this stretch are stored in passing through curcounts. We also assume that all 'ischanged' 
	 * is false, and waschanged is empty (for sake of speed, it won't bother actually checking them).
	 */
	curab=next.top().anchor;
	failed=false;

	do {
		curtb=next.top().target;
		if (curtb > lbin || curtb < fbin) { throw std::runtime_error("target bin index is out the specified range");}
		curdex=curtb-fbin;

		// Checking whether we need to set up a new row, or whether it's already in use.
		if (!ischanged[curdex]) {
			waschanged.push_back(curdex);
			ischanged[curdex]=true;
			curdex*=nlibs;
			for (lib=0; lib<nlibs; ++lib) { curcounts[curdex+lib]=0; }
		} else {
			curdex*=nlibs;
		}

		// Inner loop, to avoid multiple searches when the next batch of rows are in the same bin pair.
		do {
			curlib=next.top().library;
			int& libdex=indices[curlib];
			++(curcounts[curdex+curlib]);
			next.pop();
		 	if ((++libdex) < nums[curlib]) {
				next.push(coord(bptr[aptrs[curlib][libdex]], bptr[tptrs[curlib][libdex]], curlib));
			} else if (next.empty()) { 
				failed=true;
				break; 
			}
			if ( curab != next.top().anchor ) {
				failed=true;
				break;
			}
		} while ( curtb == next.top().target );
	} while (!failed);

	// Sorting so all targets are in ascending order during addition (anchor sorting is implicit).
	std::sort(waschanged.begin(), waschanged.end());
	return;
}

bool binner::empty() const { return next.empty(); }

int binner::get_nlibs() const { return nlibs; }

int binner::get_anchor() const { return curab; }

/*******************************************
 *********** Actual R interface ************
 *******************************************/
 
SEXP count_patch(SEXP all, SEXP bin, SEXP filter, SEXP firstbin, SEXP lastbin) try {
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int f=asInteger(filter);

	// Getting the indices of the first and last bin on the target chromosome.
	if (!isInteger(firstbin) || LENGTH(firstbin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int fbin=asInteger(firstbin);
	if (!isInteger(lastbin) || LENGTH(lastbin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
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
