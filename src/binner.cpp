#include "read_count.h"

int setup_pair_data (const Rcpp::List pairs, std::vector<Rcpp::IntegerVector>& anchor1, 
        std::vector<Rcpp::IntegerVector>& anchor2, std::vector<int>& nums, std::vector<int>& indices) {

    int nlibs=pairs.size();
    anchor1.resize(nlibs);
	anchor2.resize(nlibs);
	indices.resize(nlibs);
	nums.resize(nlibs);
	
	for (int i=0; i<nlibs; ++i) {
        const Rcpp::List current=pairs[i];
		if (current.size()!=2) { 
			throw std::runtime_error("interactions must be supplied as a data.frame with anchor.id and target.id"); 
        }

		// We assume anchor1, anchor2 have been ordered on R's side.
        for (int j=0; j<2; ++j) {
            const Rcpp::IntegerVector curvec(current[j]);
			switch (j) {
				case 0: 
					anchor1[i]=curvec;
					nums[i]=curvec.size();
					break;
				case 1: 
					anchor2[i]=curvec; 
					if (curvec.size()!=nums[i]) { 
                        throw std::runtime_error("vectors should be the same length"); 
                    }
					break;
				default: break;
			}
		}
	}
    return nlibs;
}

binner::binner(SEXP all, SEXP bin, int f, int l) : fbin(f), lbin(l), nbins(l-f+1), binid(bin), curab(-1), ischanged(nbins, false) {
    if (nbins <= 0) { throw std::runtime_error("number of bins must be positive"); }
    nlibs=setup_pair_data(all, anchor1, anchor2, nums, indices);

    // Populating the priority queue.(1-based indexing).
	for (int i=0; i<nlibs; ++i) {
        if (nums[i]) { 
            next.push(coord(binid[anchor1[i][0]-1], binid[anchor2[i][0]-1], i)); 
        }
    }

    curcounts.resize(nbins*nlibs);
	return;
}

void binner::fill() { 
    /* Resetting 'ischanged' (which indicates whether we need to set all counts to zero) and 
     * 'waschanged' (which provides a fast way to get to true values of 'ischanged').
     */
    for (const auto& wc : waschanged) {
        ischanged[wc]=false;
    }
    waschanged.clear();

	/* Running through all libraries. The idea is to use stretches of identical bin anchors, such that
	 * we only need to worry about different bin targets (i.e., the problem becomes 1-dimensional).
	 * This assumes that the bin transformation is monotonic, and that anchors are sorted. The function 
	 * will stop once all identical bin anchors have been processed. Counts for all bins 
	 * in this stretch are stored in passing through curcounts. 
	 */
	curab=next.top().anchor;
	bool failed=false;

	do {
	    const int curtb=next.top().target;
		if (curtb > lbin || curtb < fbin) { throw std::runtime_error("target bin index is out the specified range");}
		int curdex=curtb-fbin;

		// Checking whether we need to set up a new row, or whether it's already in use.
		if (!ischanged[curdex]) {
			waschanged.push_back(curdex);
			ischanged[curdex]=true;
			curdex*=nlibs;
            std::fill(curcounts.begin()+curdex, curcounts.begin()+curdex+nlibs, 0);
		} else {
			curdex*=nlibs;
		}

		// Inner loop, to avoid multiple searches when the next batch of rows are in the same bin pair.
		do {
			const int curlib=next.top().library;
			int& libdex=indices[curlib];
			++(curcounts[curdex+curlib]);
			next.pop();
		 	if ((++libdex) < nums[curlib]) {
                next.push(coord(binid[anchor1[curlib][libdex]-1], binid[anchor2[curlib][libdex]-1], curlib)); // 1-based indexing.
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

const std::vector<int>& binner::get_counts() const { return curcounts; }

const std::deque<int>& binner::get_changed()  const { return waschanged; }
