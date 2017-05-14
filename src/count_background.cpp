#include "read_count.h"
#include "neighbors.h"

SEXP count_background(SEXP all, SEXP bin, SEXP back_width, SEXP exclude_width, SEXP filter, 
		SEXP first_target_bin, SEXP last_target_bin, SEXP first_anchor_bin, SEXP last_anchor_bin) try {

	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int f=asInteger(filter);
	if (!isInteger(back_width) || LENGTH(back_width)!=1) { throw std::runtime_error("width of neighborhood regions must be an integer scalar"); }
	const int bwidth=asInteger(back_width);
	if (!isInteger(exclude_width) || LENGTH(exclude_width)!=1) { throw std::runtime_error("exclusion width must be an integer scalar"); }
	const int xwidth=asInteger(exclude_width);

	// Getting the indices of the first and last bin on the target and anchor chromosomes.
	if (!isInteger(first_target_bin) || LENGTH(first_target_bin)!=1) { throw std::runtime_error("index of first bin on target chromosome must be an integer scalar"); }
	const int ftbin=asInteger(first_target_bin);
	if (!isInteger(last_target_bin) || LENGTH(last_target_bin)!=1) { throw std::runtime_error("index of last bin on target chromosome must be an integer scalar"); }
	const int ltbin=asInteger(last_target_bin);
	if (!isInteger(first_anchor_bin) || LENGTH(first_anchor_bin)!=1) { throw std::runtime_error("index of first bin on anchor chromosome must be an integer scalar"); }
	const int fabin=asInteger(first_anchor_bin);
	if (!isInteger(last_anchor_bin) || LENGTH(last_anchor_bin)!=1) { throw std::runtime_error("index of last bin on anchor chromosome must be an integer scalar"); }
	const int labin=asInteger(last_anchor_bin);
	const bool intra=(fabin==ftbin);

	// Setting up the binning engine.
	binner engine(all, bin, ftbin, ltbin);
	const int nlibs=engine.get_nlibs();

	// Setting up the memory containers.
	const int ntbins=ltbin-ftbin+1, nabins=labin-fabin+1;
	std::deque<int> counts, anchors, targets;

	// Stuff required to compute the neighborhood.
	std::deque<int> ref_targets, ref_anchors;
	std::deque<std::deque<int> > ref_count(nlibs);
	const size_t nmodes=4, startmode=(intra ? 0 : 1);
    size_t mode;
    std::deque<std::deque<int> > neighborcount(nmodes), neighborarea(nmodes);

    basic * base_ptr=NULL;
	size_t saved_dex=0; // Points to the entry in 'anchors' to be interrogated.
	int ref_matching_dex=0; // Points to entry in 'ref_anchors' matching anchor for previous 'saved_dex'.
	size_t n_to_drop;

	// Other assorted sundries.
	int rowdex, countsum, curlib, curanchor;
	int leftbound, rightbound, leftdex, rightdex, desired_anchor;
   	size_t saved_copy_dex;
    std::deque<int>::const_iterator ccIt, wcIt;

	while (1) { 
		if (!engine.empty()) {
			engine.fill();
			curanchor=engine.get_anchor() - fabin;
            const std::deque<int>& waschanged=engine.get_changed();
            const std::deque<int>& curcounts=engine.get_counts();

            for (wcIt=waschanged.begin(); wcIt!=waschanged.end(); ++wcIt) {
				rowdex=(*wcIt)*nlibs;
                ccIt=curcounts.begin()+rowdex;

                // Adding values for neighbourhood calculations.
				ref_anchors.push_back(curanchor);
				ref_targets.push_back(*wcIt);
                countsum=0;
				for (curlib=0; curlib<nlibs; ++curlib, ++ccIt) { 
                    const int& curcount=*ccIt;
                    ref_count[curlib].push_back(curcount);
                    countsum+=curcount;
                }

			    // Adding it to the main list, if it's large enough.
				if (countsum >= f) { 
					anchors.push_back(curanchor);
					targets.push_back(ref_targets.back());
                    ccIt-=nlibs;
					for (curlib=0; curlib<nlibs; ++curlib, ++ccIt) { 
                        counts.push_back(*ccIt);
                    }
				}
			}

			// Matching sizes of the vectors.
			for (mode=startmode; mode<nmodes; ++mode) { 
				neighborarea[mode].resize(anchors.size());
				neighborcount[mode].resize(counts.size()); 
			} 
		}

        /* Computing neighbourhoods if there are enough anchors collected to define the neighbourhood of 'saved_dex',
         * or if no new information is to be added from the binning engine.
         */
		if (saved_dex < anchors.size() && (engine.empty() || ref_anchors.back() >= anchors[saved_dex] + bwidth)) { 
			bottomright br(bwidth, ntbins, intra, xwidth);
			leftright lr(bwidth, ntbins, intra, xwidth);
			updown ud(bwidth, ntbins, intra, xwidth);
			allaround aa(bwidth, ntbins, intra, xwidth);

            const int& saved_anchor=anchors[saved_dex];
			int fullsize=ref_anchors.size();
            std::deque<int> temp(nlibs), last_temp(nlibs);

			// Computing the neighborhood count for all bin pairs with anchors of 'anchors[saved_dex]'.
			for (mode=startmode; mode<nmodes; ++mode) { 
				leftdex=rightdex=0;
                std::deque<int>& curarea=neighborarea[mode];
                std::deque<int>& curcounts=neighborcount[mode];
                std::fill(temp.begin(), temp.end(), 0);

				switch (mode) {
					case 0: base_ptr=&br; break;
					case 1: base_ptr=&ud; break;
					case 2: 
                        base_ptr=&lr; 
                        leftdex=rightdex=ref_matching_dex; // starting from the stretch in ref_anchors that is past the previous 'saved_anchor'.
                        break;
					case 3: base_ptr=&aa; break;	
				}

				/* If a bump_level results in an increase in desired_anchor, we can hot-start 
				 * from the existing left/rightdex from the last level (i.e., no need to reset to 
				 * zero within the loop). Both indices MUST increase (even when saved_copy_dex 
				 * is reset to saved_dex) as they've been stuck on the previous desired_anchor.
                 *
                 * If bumping doesn't result in an increase, we need to hot-start from the 
                 * previous left/rightdex that we used for the current level. This necessitates
                 * storing of the last_anchor and the last left/right indices. We also need 
                 * to store the temporary counts corresponding to those indices.
                 *
                 * Note that we reset to zero in quadrant_bg, which needs to calculate the 
                 * neighbourhood across _all_ pairs; not just those at a particular anchor.
				 */
                int last_anchor=-1, last_leftdex=0, last_rightdex=0;
				do { 
   	        		base_ptr->set(saved_anchor, 0); // Dummy input to get desired_anchor.
                    desired_anchor=base_ptr->row;
				    if (desired_anchor < 0 || desired_anchor >= nabins) { continue; }
                    if (last_anchor==desired_anchor) {
                        leftdex=last_leftdex;
                        rightdex=last_rightdex;
                        std::copy(last_temp.begin(), last_temp.end(), temp.begin());
                    } else {
                        last_anchor=desired_anchor;
                        last_leftdex=leftdex;
                        last_rightdex=rightdex;
                        std::copy(temp.begin(), temp.end(), last_temp.begin());
                    }
                    
					for (saved_copy_dex=saved_dex; saved_copy_dex < anchors.size() && 
							anchors[saved_copy_dex]==saved_anchor; ++saved_copy_dex) { 
						base_ptr->set(saved_anchor, targets[saved_copy_dex]);
						leftbound=base_ptr->left;
						rightbound=base_ptr->right;

						// Identifying all reference elements associated with the neighborhood of this saved_copy_dex at the current level.
						while (rightdex < fullsize && (ref_anchors[rightdex] < desired_anchor || 
									(ref_anchors[rightdex]==desired_anchor && ref_targets[rightdex] < rightbound))) { 
                            for (curlib=0; curlib<nlibs; ++curlib) { temp[curlib]+=ref_count[curlib][rightdex]; }
							++rightdex;
						}
						while (leftdex < fullsize && (ref_anchors[leftdex] < desired_anchor || 
									(ref_anchors[leftdex]==desired_anchor && ref_targets[leftdex] < leftbound))) {
                            for (curlib=0; curlib<nlibs; ++curlib) { temp[curlib]-=ref_count[curlib][leftdex]; }
							++leftdex;
						}
		
						curarea[saved_copy_dex]+=rightbound-leftbound;
                        rowdex=saved_copy_dex*nlibs;
                        for (curlib=0; curlib<nlibs; ++curlib) { curcounts[rowdex+curlib]+=temp[curlib]; }
					}
				} while (base_ptr->bump_level());

				/* Storing the location on ref_anchors where 'saved_anchor' terminates, to hotstart for
				 * the next 'saved_anchor' for 'leftright' (as it only operates within the same anchor values).
				 */
				if (mode==2) { ref_matching_dex=rightdex; }
			}
			
            // Shifting onwards to the next 'saved_anchor'.
			while (saved_dex < anchors.size() && anchors[saved_dex]==saved_anchor) { ++saved_dex; } 
		}

		// Dropping elements in the reference data that will no longer be used, to save memory.
		if (saved_dex < anchors.size()) { 
			n_to_drop=0;
			while (n_to_drop < ref_anchors.size() && ref_anchors[n_to_drop] < anchors[saved_dex] - bwidth) { ++n_to_drop; }
			if (n_to_drop) { 
				ref_anchors.erase(ref_anchors.begin(), ref_anchors.begin()+n_to_drop);
				ref_targets.erase(ref_targets.begin(), ref_targets.begin()+n_to_drop);
                for (curlib=0; curlib<nlibs; ++curlib) { 
                    std::deque<int>& cur_count=ref_count[curlib];
                    cur_count.erase(cur_count.begin(), cur_count.begin()+n_to_drop);
                }

				/* The 'matching_dex' gets pulled down as well. It's possible that the next 'saved_anchor'
				 * skips a couple of ref_anchors, so the number dropped might include 'matching_dex'. In
				 * that case, we just start from zero during the neighborhood calculations.
				 */
				ref_matching_dex-=n_to_drop;
				if (ref_matching_dex < 0) { ref_matching_dex=0; }
			}
		} else if (engine.empty()) { break; }
	}

	// Storing into R objects.
	SEXP output=PROTECT(allocVector(VECSXP, 5));
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
	
        // Storing the count statistics for each neighbourhood.
        SET_VECTOR_ELT(output, 3, allocVector(VECSXP, nmodes));
        SEXP count_out=VECTOR_ELT(output, 3); 
        std::deque<std::deque<int*> > noptrs(nmodes);
        for (mode=0; mode<nmodes; ++mode) {
            SET_VECTOR_ELT(count_out, mode, allocMatrix(INTSXP, ncombos, nlibs));
            noptrs[mode].resize(nlibs);
            for (curlib=0; curlib<nlibs; ++curlib) {
                if (curlib==0) { noptrs[mode][0]=INTEGER(VECTOR_ELT(count_out, mode)); }
                else { noptrs[mode][curlib]=noptrs[mode][curlib-1] + ncombos; }
            }
        }

        // Setting up output for the neighbourhood size.
        SET_VECTOR_ELT(output, 4, allocVector(VECSXP, nmodes));
        SEXP n_out=VECTOR_ELT(output, 4); 
        std::deque<int*> nnptrs(nmodes);
        for (mode=0; mode<nmodes; ++mode) {
            SET_VECTOR_ELT(n_out, mode, allocVector(INTSXP, ncombos));
            nnptrs[mode]=INTEGER(VECTOR_ELT(n_out, mode));
        }
		
        // Filling in empty info for inter-chromosomal interactions.
        if (startmode!=0) {
            for (curlib=0; curlib<nlibs; ++curlib) {
                std::fill(noptrs[0][curlib], noptrs[0][curlib]+ncombos, 0);
            }
            std::fill(nnptrs[0], nnptrs[0]+ncombos, 0);
        }

		// Iterating across and filling both the matrix and the components.
		int cdex=-1;
		double backaverage, maxaverage;
		for (size_t vecdex=0; vecdex<anchors.size(); ++vecdex) {
			aoptr[vecdex]=anchors[vecdex]+fabin;
			toptr[vecdex]=targets[vecdex]+ftbin;

            for (curlib=0; curlib<nlibs; ++curlib) {
				++cdex; 			   	
				coptrs[curlib][vecdex]=counts[cdex]; 

                for (mode=startmode; mode<nmodes; ++mode) {
                    noptrs[mode][curlib][vecdex]=neighborcount[mode][cdex];
                }
			}
        }

        for (mode=startmode; mode<nmodes; ++mode) {
            std::copy(neighborarea[mode].begin(), neighborarea[mode].end(), nnptrs[mode]);
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


