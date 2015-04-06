#include "diffhic.h"
#include "coord.h"

SEXP count_patch(SEXP all, SEXP bin, SEXP filter) try {
    if (!isInteger(bin)) { throw std::runtime_error("anchor bin indices must be integer vectors"); }
	const int* bptr=INTEGER(bin)-1; // Assuming 1-based indices for anchors and targets.
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int f=asInteger(filter);
   	if (!isNewList(all)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
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
            int* ptr=INTEGER(current_col);
            switch (j) {
                case 0: 
					aptrs[i]=ptr; 
					nums[i]=LENGTH(current_col);
					break;
                case 1: tptrs[i]=ptr; break;
                default: break;
            }
		}
		
		// Populating the priority queue.
		if (nums[i]) { next.push(coord(bptr[aptrs[i][0]], bptr[tptrs[i][0]], i)); }
	}
	
	/* Running through all libraries. The idea is to use stretches of identical anchors to restrain the map
	 * so that it only contains differing targets. This assumes that the bin transformation is monotonic,
	 * and that anchors are sorted.
	 */
	std::deque<int> counts, anchors, targets;
	int total=0, countsum=0;
	std::map<int, int> bins;
	std::map<int, int>::iterator itb;
	std::deque<int> curcounts;

	int curab, curtb, curlib;
	bool failed;
	while (!next.empty()) {
		curab=next.top().anchor;
		failed=false;
		do {
			curtb=next.top().target;
			itb=bins.lower_bound(curtb);
			if (itb==bins.end() || bins.key_comp()(curtb, itb->first)) {
				curcounts.resize(total+nlibs);
				itb=bins.insert(itb, std::make_pair(curtb, total));
				total+=nlibs;
			}
			const int& curdex=itb->second;

			// Inner loop, to avoid multiple searches when the next batch of rows are in the same bin.
			do {
				curlib=next.top().library;
				int& libdex=indices[curlib];
				curcounts[curdex+curlib]+=1;
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

		// Adding it to the main list, if it's large enough.
		for (itb=bins.begin(); itb!=bins.end(); ++itb) {
			const int& curdex=itb->second;
			countsum=0;
			for (int i=0; i<nlibs; ++i) { countsum+=curcounts[curdex+i]; }
			if (countsum >= f) { 
				anchors.push_back(curab);
				targets.push_back(itb->first);
				for (int i=0; i<nlibs; ++i) { counts.push_back(curcounts[curdex+i]); }
			}	
		}
		total=0;
		bins.clear();
		curcounts.clear();
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
