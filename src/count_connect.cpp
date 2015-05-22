#include "read_count.h"

SEXP count_connect(SEXP all, SEXP start, SEXP end, SEXP region, SEXP filter, SEXP is_second) try {		
	if (!isInteger(start) || !isInteger(end) || !isInteger(region)) { 
		throw std::runtime_error("fragment/region indices must be integer vectors"); }
	const int ni=LENGTH(start), nrp1=LENGTH(region)+1; // As 'end' refers to one-past-the-end. 
	if (LENGTH(end)!=ni) { throw std::runtime_error("start/end index vectors should be the same length"); }

	// Setting up pointers (-1 for 1-based indexing).
	const int * sptr=INTEGER(start)-1,
			* eptr=INTEGER(end)-1,
			* rptr=INTEGER(region)-1;

	// Checking scalars.
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int filtval=asInteger(filter);

	// Checking second specifier.
	const bool check_second=(is_second!=R_NilValue);
	const int* iptr=NULL;
	if (check_second) {
 	   if (!isLogical(is_second)) { throw std::runtime_error("secondary specifier should be a logical vector or NULL"); }
	   iptr=LOGICAL(is_second)-1;
	   const int nseconds=LENGTH(is_second);
	   for (int it=1; it<nrp1; ++it) {
		   if (rptr[it] <= 0 || rptr[it] > nseconds) { throw std::runtime_error("region indices out of range for secondary specifier vector"); }
	   }
	}

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
   	if (!isNewList(all)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }
	const int nlibs=LENGTH(all);
	std::deque<const int*> aptrs(nlibs), tptrs(nlibs);
	std::deque<int> nums(nlibs), indices(nlibs);
	std::priority_queue<coord, std::deque<coord>, std::greater<coord> > next;

	for (int i=0; i<nlibs; ++i) {
		SEXP current=VECTOR_ELT(all, i);
		if (!isNewList(current) || LENGTH(current)!=2) { 
			throw std::runtime_error("interactions must be supplied as a data.frame with anchor.id, target.id"); }

		for (int j=0; j<2; ++j) {
			SEXP current_col=VECTOR_ELT(current, j);
			if (!isInteger(current_col)) { throw std::runtime_error("interaction data must be in integer format"); }
			int* ptr=INTEGER(current_col);
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
		if (nums[i]) { next.push(coord(aptrs[i][0], tptrs[i][0], i)); }
	}
	
	// Running through all libraries.
	std::deque<int> counts;
	typedef std::pair<int, int> combo;
	std::map<combo, std::pair<int, int> > bins;
	std::map<combo, std::pair<int, int> >::iterator itb;
	std::deque<int> curcounts(nlibs);

	int curab, curtb, curlib;
	combo temp;
	int counter=0;
	int x1, x2, lib;

	int smallest=-1;
	std::deque<int> returned_anchors, returned_targets;
	std::deque<size_t> count_indices;
	int countsum;
	size_t index;

	while (!next.empty()) {
		curab=next.top().anchor;
		curtb=next.top().target;
		++counter;

		// This speeds things up by collecting entries with the same fragment indices.
		do {
			curlib=next.top().library;
			int& libdex=indices[curlib];
			curcounts[curlib]+=1;
			next.pop();
			if ((++libdex) < nums[curlib]) {
				next.push(coord(aptrs[curlib][libdex], tptrs[curlib][libdex], curlib));
			} 
		} while (!next.empty() && next.top().anchor==curab && next.top().target==curtb);

		/* Allocating counts to every pair of ranges containing these fragments. This
 		 * shouldn't be too sadistic if there aren't too many overlapping ranges.
 		 */
		if (curab>ni) { throw std::runtime_error("invalid anchor index for supplied fragments"); } // 1-based indexing, so '>' is right.
		const int& s1x=sptr[curab];
		const int& e1x=eptr[curab];
		if (curtb>ni) { throw std::runtime_error("invalid target index for supplied fragments"); }
		const int& s2x=sptr[curtb];
		const int& e2x=eptr[curtb];
	
		if (s1x!=e1x && s2x!=e2x) { 
			if (s1x <= 0 || s2x <= 0 || e1x > nrp1 || e2x > nrp1) { throw std::runtime_error("invalid start/endpoints for region indices"); }
			for (x1=s1x; x1<e1x; ++x1) {
				for (x2=s2x; x2<e2x; ++x2) { 
					if (rptr[x1] > rptr[x2]) {
						temp.first=rptr[x1];
						temp.second=rptr[x2];
					} else {
						// Avoid redundant naming, when looking within the same ranges.
						temp.first=rptr[x2];
						temp.second=rptr[x1];
					}
					if (temp.first < smallest) { throw std::runtime_error("largest index cannot be lower than 'smallest' threshold"); } 

					// Checking secondary; i.e., one read in second ranges, other region outside.
					if (check_second && iptr[temp.first]==iptr[temp.second]) { continue; }

					itb=bins.lower_bound(temp);
					if (itb==bins.end() || bins.key_comp()(temp, itb->first)) {
						itb = bins.insert(itb, std::make_pair(temp, std::make_pair(counts.size(), counter)));
						counts.resize(counts.size()+nlibs);
					} else if ((itb->second).second==counter) { 
						/* The 'counter' avoids adding the same range twice to a particular pair. Regions
						 * can be irregularly sized and spaced, e.g., nested, so it's not possible to set up 
						 * a general rule to avoid redundant counting. For example, the target fragment might
						 * overlap regions 1, 2 while the anchor might just overlap 1, if 2 is nested in 1.
						 */
						continue;
					}
					(itb->second).second=counter;
					const int& index=(itb->second).first;
					for (lib=0; lib<nlibs; ++lib) { counts[index+lib] += curcounts[lib]; }
				}
			}
		}

		// Resetting.
		for (lib=0; lib<nlibs; ++lib) { curcounts[lib]=0; }

		/* Removing all entries where the first index is below the smallest region index for the anchor 
		 * fragment. As we iterate, the anchor fragment will increase, so the smallest region index cannot 
		 * decrease (assuming regions have already been sorted at the R-level). This means that lower entries 
		 * will no longer be accessible. Swapping with the index for the target fragment has no effect, as 
		 * you need to be higher to swap. The aim is to shrink the map to make it (almost) 1-dimensional.
		 */
		smallest=rptr[s1x];
		for (x1=s1x+1; x1<e1x; ++x1) { if (smallest > rptr[x1]) { smallest=rptr[x1]; } }
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
		for (int i=0; i<nlibs; ++i) {
			if (i==0) { coptrs[i]=INTEGER(VECTOR_ELT(output, 2)); }
			else { coptrs[i]=coptrs[i-1]+ncombos; }
		}	
		
		// Iterating across and filling both the matrix and the components.
		for (size_t odex=0; odex < ncombos; ++odex) { 
			aoptr[odex]=returned_anchors[odex];
			toptr[odex]=returned_targets[odex];
			const int& index=count_indices[odex];
			for (int i=0; i<nlibs; ++i) { coptrs[i][odex]=counts[index+i]; }
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
