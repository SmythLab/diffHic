#include "diffhic.h"
#include "coord.h"

SEXP count_patch(SEXP all, SEXP bin, SEXP filter, SEXP firstbin, SEXP lastbin) try {
	if (!isInteger(bin)) { throw std::runtime_error("anchor bin indices must be integer vectors"); }
	const int* bptr=INTEGER(bin)-1; // Assuming 1-based indices for anchors and targets.
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int f=asInteger(filter);
	if (!isInteger(firstbin) || LENGTH(firstbin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }

	// Getting the indices of the first and last bin on the target chromosome.
	const int fbin=asInteger(firstbin);
	if (!isInteger(lastbin) || LENGTH(lastbin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int lbin=asInteger(lastbin);

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
	if (!isNewList(all)) { throw std::runtime_error("data on interacting read pairs must be contained within a list"); }
	const int nlibs=LENGTH(all);
	std::deque<const int*> aptrs(nlibs), tptrs(nlibs);
	std::deque<int> nums(nlibs), indices(nlibs);
	std::priority_queue<coord, std::deque<coord>, std::greater<coord> > next;
	
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

	// Setting up the memory containers.
	const int nbins=lbin-fbin+1;
	int* curcounts=(int*)R_alloc(nlibs*nbins, sizeof(int)); 
	bool* ischanged=(bool*)R_alloc(nbins, sizeof(bool));
	for (int i=0; i<nbins; ++i) { ischanged[i]=false; }
	std::deque<int> waschanged, counts, anchors, targets;
	size_t rowdex;

	/* Running through all libraries. The idea is to use stretches of identical bin anchors, such that
	 * we only need to worry about different bin targets (i.e., the problem becomes 1-dimensional).
	 * This assumes that the bin transformation is monotonic, and that anchors are sorted.
	 */
	int countsum=0;
	int curab, curtb, curlib, curdex, lib;
	bool failed;

	while (!next.empty()) {
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

		// Adding it to the main list, if it's large enough.
		for (rowdex=0; rowdex<waschanged.size(); ++rowdex) {
			curdex=waschanged[rowdex]*nlibs;
			countsum=0;
			for (lib=0; lib<nlibs; ++lib) { countsum+=curcounts[curdex+lib]; }
			if (countsum >= f) { 
				anchors.push_back(curab);
				targets.push_back(waschanged[rowdex]+fbin);
				for (lib=0; lib<nlibs; ++lib) { counts.push_back(curcounts[curdex+lib]); }
			}

			// Resetting the ischanged vector for the next stretch of bin anchors.
			ischanged[waschanged[rowdex]]=false;
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
		for (int i=0; i<nlibs; ++i) {
			if (i==0) { coptrs[i]=INTEGER(VECTOR_ELT(output, 2)); }
			else { coptrs[i]=coptrs[i-1]+ncombos; }
		}	
		
		// Iterating across and filling both the matrix and the components.
		int cdex=-1;
		for (size_t i=0; i<anchors.size(); ++i) {
			aoptr[i]=anchors[i];
			toptr[i]=targets[i];
			for (int j=0; j<nlibs; ++j) { coptrs[j][i]=counts[++cdex]; }
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
