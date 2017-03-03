#include "read_count.h"

binner::binner(SEXP all, SEXP bin, int f, int l) : fbin(f), lbin(l), nbins(l-f+1), ischanged(NULL), curcounts(NULL) {
	if (!isInteger(bin)) { throw std::runtime_error("anchor bin indices must be integer vectors"); }
	bptr=INTEGER(bin)-1; // Assuming 1-based indices for anchors and targets.
    if (nbins <= 0) { throw std::runtime_error("number of bins must be positive"); }

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

    try {
        ischanged=new bool[nbins];
        std::fill(ischanged, ischanged+nbins, false);
        curcounts=new int[nbins*nlibs];
    } catch (std::exception& e) {
        delete [] ischanged;
        delete [] curcounts;
        throw;
    }
	return;
}

binner::~binner () {
    delete [] ischanged;
    delete [] curcounts;
    return;
}

void binner::fill() { 
    /* Resetting 'ischanged' (which indicates whether we need to set all counts to zero) and 
     * 'waschanged' (which provides a fast way to get to true values of 'ischanged').
     */
    for (changedex=0; changedex<waschanged.size(); ++changedex) {
        ischanged[waschanged[changedex]]=false;
    }
    waschanged.clear();

	/* Running through all libraries. The idea is to use stretches of identical bin anchors, such that
	 * we only need to worry about different bin targets (i.e., the problem becomes 1-dimensional).
	 * This assumes that the bin transformation is monotonic, and that anchors are sorted. The function 
	 * will stop once all identical bin anchors have been processed. Counts for all bins 
	 * in this stretch are stored in passing through curcounts. 
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
            std::fill(curcounts+curdex, curcounts+curdex+nlibs, 0);
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

int binner::get_nbins() const { return nbins; }

int binner::get_anchor() const { return curab; }

const int* binner::get_counts() const { return curcounts; }

const std::deque<int>& binner::get_changed()  const { return waschanged; }
