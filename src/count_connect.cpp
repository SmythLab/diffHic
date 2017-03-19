#include "read_count.h"

void accumulate_pairs(const int& curab, const int& curtb, pair_queue& next, 
        std::deque<const int*>& aptrs, std::deque<const int*>& tptrs, 
        std::deque<int>& nums, std::deque<int>& indices, 
        std::deque<int>& counts) {
    // This speeds things up by collecting entries with the same fragment indices.
    int curlib;
    do {
        curlib=next.top().library;
        int& libdex=indices[curlib];
        counts[curlib]+=1;
        next.pop();
        if ((++libdex) < nums[curlib]) {
            next.push(coord(aptrs[curlib][libdex], tptrs[curlib][libdex], curlib));
        } 
    } while (!next.empty() && next.top().anchor==curab && next.top().target==curtb);
    return;
}

SEXP count_connect(SEXP all, SEXP start1, SEXP end1, SEXP region1, 
        SEXP start2, SEXP end2, SEXP region2, SEXP filter) try {		
   
	if (!isInteger(start1) || !isInteger(end1)) { throw std::runtime_error("fragment indices (1) must be integer vectors"); }
	const int ni1=LENGTH(start1);
	if (LENGTH(end1)!=ni1) { throw std::runtime_error("start/end index vectors (1) should be the same length"); }

    if (!isInteger(region1)) { throw std::runtime_error("region indices (1) must be integer vectors"); }
    const int nrp1=LENGTH(region1)+1; // As 'end' refers to one-past-the-end. 

    // Setting up pointers (-1 for 1-based indexing).
	const int * s1ptr=INTEGER(start1)-1,
			* e1ptr=INTEGER(end1)-1,
			* r1ptr=INTEGER(region1)-1,
            * s2ptr=s1ptr, * e2ptr=e1ptr, * r2ptr=r1ptr;
    
    // Checking scalars.
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int filtval=asInteger(filter);

    // Repeating for the second region, if specified.
    int ni2=ni1, nrp2=nrp1;
    const bool use_second=(region2!=R_NilValue);
    if (use_second) {
        if (!isInteger(start2) || !isInteger(end2)) {
            throw std::runtime_error("fragment indices (2) must be integer vectors"); }
        ni2=LENGTH(start2); 
        if (LENGTH(end2)!=ni2) { throw std::runtime_error("start/end index vectors (2) should be the same length"); }

        if (!isInteger(region2)) { throw std::runtime_error("region indices (2) must be integer vectors"); }
        nrp2=LENGTH(region2)+1; 

        s2ptr=INTEGER(start2)-1;
        e2ptr=INTEGER(end2)-1;
        r2ptr=INTEGER(region2)-1;
    }

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
   	if (!isNewList(all)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }
	const int nlibs=LENGTH(all);
	std::deque<const int*> aptrs(nlibs), tptrs(nlibs);
	std::deque<int> nums(nlibs), indices(nlibs);
    setup_pair_data(all, nlibs, aptrs, tptrs, nums, indices);
    
    pair_queue next;
    int lib;
    for (lib=0; lib<nlibs; ++lib) {
        if (nums[lib]) { next.push(coord(aptrs[lib][0], tptrs[lib][0], lib)); }
    }
	
	// Running through all libraries.
	std::deque<int> counts;
	typedef std::pair<int, int> combo;
	std::map<combo, std::pair<int, int> > bins;
	std::map<combo, std::pair<int, int> >::iterator itb;
	std::deque<int> curcounts(nlibs);

	int curab, curtb;
	combo temp;
	int mode, counter=0;
	int x1, x2, y1, y2;

	int smallest=-1;
	std::deque<int> returned_anchors, returned_targets;
	std::deque<size_t> count_indices;
	int countsum;
	size_t index;

	while (!next.empty()) {
		curab=next.top().anchor;
		curtb=next.top().target;
		++counter;
        accumulate_pairs(curab, curtb, next, aptrs, tptrs, nums, indices, curcounts); 

		/* Allocating counts to every pair of ranges containing these fragments. This
         * has a lot of loops but it shouldn't be too bad if there aren't many 
         * overlapping ranges. The 'uppermode' just checks the anchor vs the second
         * set of regions and the target vs the first set of regions, if there are two sets.
 		 */
        for (mode=0; mode<(use_second ? 2 : 1); ++mode) { 
            if (mode==0) { 
                y1=curab;
                y2=curtb;
            } else {
                y1=curtb;
                y2=curab;
            }
            if (y1>ni1) { throw std::runtime_error("invalid index (1) for supplied fragments"); } // 1-based indexing, so '>' is right.
            const int& s1x=s1ptr[y1];
            const int& e1x=e1ptr[y1];
            if (y2>ni2) { throw std::runtime_error("invalid index (2) for supplied fragments"); }
            const int& s2x=s2ptr[y2];
            const int& e2x=e2ptr[y2];
        
            if (s1x!=e1x && s2x!=e2x) { 
                if (s1x <= 0 || s2x <= 0 || e1x > nrp1 || e2x > nrp2) { throw std::runtime_error("invalid start/endpoints for region indices"); }
                for (x1=s1x; x1<e1x; ++x1) {
                    for (x2=s2x; x2<e2x; ++x2) {

                        /* Avoid redundant naming, when looking within the same ranges; region with 
                         * higher genomic coordinates goes first. This should always increase, see below.
                         */
                        if (r1ptr[x1] > r2ptr[x2]) {
                            temp.first=r1ptr[x1];
                            temp.second=r2ptr[x2];
                        } else {
    						temp.first=r2ptr[x2];
    						temp.second=r1ptr[x1];
    					}
	    				if (temp.first < smallest) { throw std::runtime_error("largest index cannot be lower than 'smallest' threshold"); }

    					itb=bins.lower_bound(temp);
    					if (itb==bins.end() || bins.key_comp()(temp, itb->first)) {
    						itb = bins.insert(itb, std::make_pair(temp, std::make_pair(counts.size(), counter)));
    						counts.resize(counts.size()+nlibs);
    					} else if ((itb->second).second==counter) { 
    						/* The 'counter' avoids adding the same range twice to a particular pair. Regions
    						 * can be irregularly sized and spaced, e.g., nested, so it's not possible to set up 
    						 * a loop condition involving x1 and x2 to avoid redundant counting. For example, the 
                             * target fragment might overlap regions 1, 2 while the anchor might just overlap 1, 
                             * if 2 is nested in 1; so mandating that 'x1 >= x2' won't work. In addition, you 
                             * need to coordinate between the 'mode' swap of y1/2.
    						 */
    						continue;
    					}
    					(itb->second).second=counter;
    					const int& index=(itb->second).first;
    					for (lib=0; lib<nlibs; ++lib) { counts[index+lib] += curcounts[lib]; }
                    }
				}
			}
		}

		// Resetting counts.
        std::fill(curcounts.begin(), curcounts.end(), 0);

		/* Saving all entries where the first region index is below the smallest anchor-overlapped region index for the current read pair.
         * This frees up the map to make it almost 1-dimensional, which should reduce memory usage and improve speed.
         *
         * The logic is that, as the anchor index increases, the anchor fragment cannot overlap a region with a lower index in 'r(1/2)ptr'
         * (assuming that the regions have been sorted at the R level). So, the anchor-overlapped region index cannot decrease.
         * Now, the first index of each entry is the larger of the two indices; if an entry has a first index below the current
         * smallest anchor-overlapped region, it cannot possibly be updated at later iterations. This means we can remove it.
  		 */
        const int& s1x=s1ptr[curab];
        const int& e1x=e1ptr[curab];
        for (x1=s1x; x1<e1x; ++x1) { if (smallest > r1ptr[x1]) { smallest=r1ptr[x1]; } }
        if (use_second) { 
            const int& s2x=s2ptr[curab];
            const int& e2x=e2ptr[curab];
            for (x2=s2x; x2<e2x; ++x2) { if (smallest > r2ptr[x2]) { smallest=r2ptr[x2]; } }
        }

		itb=bins.begin();
		while (itb!=bins.end()) { 
			if ((itb->first).first < smallest) {
				countsum=0; 
				index=(itb->second).first;
				for (lib=0; lib<nlibs; ++lib, ++index) { countsum += counts[index]; }

				if (countsum >= filtval) { // Checking that the count sum is sufficient.
					returned_anchors.push_back((itb->first).first);
					returned_targets.push_back((itb->first).second);
					count_indices.push_back((itb->second).first);
				}
				bins.erase(itb++);
			} else {
				break;
			}
		}
	}

	// Assessing how many combinations are above threshold.
	for (itb=bins.begin(); itb!=bins.end(); ++itb) { 
		countsum=0;
		index=(itb->second).first;
		for (lib=0; lib<nlibs; ++lib, ++index) { countsum += counts[index]; }
		if (countsum >= filtval) { 
			returned_anchors.push_back((itb->first).first);
			returned_targets.push_back((itb->first).second);
			count_indices.push_back((itb->second).first);
		}
	}
	const size_t ncombos=count_indices.size();

	// Returning all count combinations underneath the threshold.
	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ncombos));
		int* aoptr=INTEGER(VECTOR_ELT(output, 0));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ncombos));
		int* toptr=INTEGER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, allocMatrix(INTSXP, ncombos, nlibs));
		std::deque<int*> coptrs(nlibs);
		for (lib=0; lib<nlibs; ++lib) {
			if (lib==0) { coptrs[lib]=INTEGER(VECTOR_ELT(output, 2)); }
			else { coptrs[lib]=coptrs[lib-1]+ncombos; }
		}	
		
		// Iterating across and filling both the matrix and the components.
		for (size_t odex=0; odex < ncombos; ++odex) { 
			aoptr[odex]=returned_anchors[odex];
			toptr[odex]=returned_targets[odex];
			const int& index=count_indices[odex];
			for (lib=0; lib<nlibs; ++lib) { coptrs[lib][odex]=counts[index+lib]; }
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

/* This does the same for DNase-C data, where the assignment into ranges 
 * has already been done at the R-level using linkOverlaps. So we just
 * need to scan across the assignments in 'links' and aggregate counts.
 */

SEXP count_reconnect(SEXP links, SEXP filter) try {	
	// Setting up the structure to the links.
 	const int nlibs=LENGTH(links);
	std::deque<const int*> aptrs(nlibs), tptrs(nlibs);
	std::deque<int> nums(nlibs), indices(nlibs);
    setup_pair_data(links, nlibs, aptrs, tptrs, nums, indices);

    pair_queue next;
    int lib;
    for (lib=0; lib<nlibs; ++lib) {
        if (nums[lib]) { next.push(coord(aptrs[lib][0], tptrs[lib][0], lib)); }
    }

    // Checking scalars.
 	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int filtval=asInteger(filter);

   // Running through the pairs.
    std::deque<int> returned_anchors, returned_targets;
	std::deque<int> counts;
	std::deque<int> curcounts(nlibs);
    int total, curab, curtb;

    while (!next.empty()) {
		curab=next.top().anchor;
		curtb=next.top().target;
        accumulate_pairs(curab, curtb, next, aptrs, tptrs, nums, indices, curcounts); 

        // Inserting if the sum of counts exceeds the value.
        total=0;
        for (lib=0; lib<nlibs; ++lib) { total += curcounts[lib]; }
        if (total >= filtval) { 
            returned_anchors.push_back(curab);
            returned_targets.push_back(curtb);
            for (lib=0; lib<nlibs; ++lib) { counts.push_back(curcounts[lib]); }
        }
        std::fill(curcounts.begin(), curcounts.end(), 0);
    }

    // Saving the output.
	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
        const int ncombos = returned_anchors.size();
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ncombos));
		int* aoptr=INTEGER(VECTOR_ELT(output, 0));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ncombos));
		int* toptr=INTEGER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, allocMatrix(INTSXP, ncombos, nlibs));
		std::deque<int*> coptrs(nlibs);
		for (lib=0; lib<nlibs; ++lib) {
			if (lib==0) { coptrs[lib]=INTEGER(VECTOR_ELT(output, 2)); }
			else { coptrs[lib]=coptrs[lib-1]+ncombos; }
		}	
		
		// Iterating across and filling both the matrix and the components.
        int index;
		for (size_t odex=0; odex < ncombos; ++odex) { 
			aoptr[odex]=returned_anchors[odex];
			toptr[odex]=returned_targets[odex];
			index=odex*nlibs;
			for (lib=0; lib<nlibs; ++lib) { coptrs[lib][odex]=counts[index+lib]; }
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

