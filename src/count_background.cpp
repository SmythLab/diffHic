#include "diffhic.h"
#include "read_count.h"
#include "neighbors.h"
#include "utils.h"

SEXP count_background(SEXP all, SEXP bin, SEXP back_width, SEXP exclude_width, SEXP filter, 
		SEXP first_anchor2_bin, SEXP last_anchor2_bin, SEXP first_anchor1_bin, SEXP last_anchor1_bin) {
    BEGIN_RCPP

    const int f=check_integer_scalar(filter, "filter value");
    const int bwidth=check_integer_scalar(back_width, "width of neighborhood regions");
    const int xwidth=check_integer_scalar(exclude_width, "exclusion width");

	// Getting the indices of the first and last bin on the anchor2 and anchor chromosomes.
    const int ftbin=check_integer_scalar(first_anchor2_bin, "index of first bin on second anchor chromosome");
    const int ltbin=check_integer_scalar(last_anchor2_bin, "index of last bin on second anchor chromosome");
    const int fabin=check_integer_scalar(first_anchor1_bin, "index of first bin on first anchor chromosome");
    const int labin=check_integer_scalar(last_anchor1_bin, "index of last bin on first anchor chromosome");
    const bool intra=(fabin==ftbin);

	// Setting up the binning engine.
	binner engine(all, bin, ftbin, ltbin);
	const int nlibs=engine.get_nlibs();

	// Setting up the memory containers.
	const int ntbins=ltbin-ftbin+1, nabins=labin-fabin+1;
	std::deque<int> counts, anchors, targets;

	// Stuff required to compute the neighborhood.
	std::deque<int> ref_targets, ref_anchors;
	std::vector<std::deque<int> > ref_count(nlibs);
	const size_t nmodes=4, startmode=(intra ? 0 : 1);
    std::vector<std::deque<int> > neighborcount(nmodes), neighborarea(nmodes);

	size_t saved_dex=0; // Points to the entry in 'anchors' to be interrogated.
	int ref_matching_dex=0; // Points to entry in 'ref_anchors' matching anchor for previous 'saved_dex'.

	while (1) { 
		if (!engine.empty()) {
			engine.fill();
			const int curanchor=engine.get_anchor() - fabin;
            const std::deque<int>& waschanged=engine.get_changed();
            const std::vector<int>& curcounts=engine.get_counts();

            for (const auto& wc : waschanged) { 
				ref_anchors.push_back(curanchor);
				ref_targets.push_back(wc);

                // Adding values for neighbourhood calculations.
				int rowdex=wc*nlibs;
                auto ccIt=curcounts.begin()+rowdex;
                int countsum=0;
				for (int curlib=0; curlib<nlibs; ++curlib, ++ccIt) { 
                    const int& curcount=*ccIt;
                    ref_count[curlib].push_back(curcount);
                    countsum+=curcount;
                }

			    // Adding it to the main list, if it's large enough.
				if (countsum >= f) { 
					anchors.push_back(curanchor);
					targets.push_back(ref_targets.back());
                    ccIt-=nlibs;
					for (int curlib=0; curlib<nlibs; ++curlib, ++ccIt) { 
                        counts.push_back(*ccIt);
                    }
				}
			}

			// Matching sizes of the vectors.
			for (size_t mode=startmode; mode<nmodes; ++mode) { 
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
			for (size_t mode=startmode; mode<nmodes; ++mode) { 
				int leftdex=0, rightdex=0;
                std::deque<int>& curarea=neighborarea[mode];
                std::deque<int>& curcounts=neighborcount[mode];
                std::fill(temp.begin(), temp.end(), 0);

                basic * base_ptr=NULL;
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
                    const int desired_anchor=base_ptr->row;
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
                    
					for (size_t saved_copy_dex=saved_dex; saved_copy_dex < anchors.size() && 
							anchors[saved_copy_dex]==saved_anchor; ++saved_copy_dex) { 
						base_ptr->set(saved_anchor, targets[saved_copy_dex]);
						const int leftbound=base_ptr->left;
						const int rightbound=base_ptr->right;

						// Identifying all reference elements associated with the neighborhood of this saved_copy_dex at the current level.
						while (rightdex < fullsize && (ref_anchors[rightdex] < desired_anchor || 
									(ref_anchors[rightdex]==desired_anchor && ref_targets[rightdex] < rightbound))) { 
                            for (int curlib=0; curlib<nlibs; ++curlib) { temp[curlib]+=ref_count[curlib][rightdex]; }
							++rightdex;
						}
						while (leftdex < fullsize && (ref_anchors[leftdex] < desired_anchor || 
									(ref_anchors[leftdex]==desired_anchor && ref_targets[leftdex] < leftbound))) {
                            for (int curlib=0; curlib<nlibs; ++curlib) { temp[curlib]-=ref_count[curlib][leftdex]; }
							++leftdex;
						}
		
						curarea[saved_copy_dex]+=rightbound-leftbound;
                        int rowdex=saved_copy_dex*nlibs;
                        for (int curlib=0; curlib<nlibs; ++curlib) { curcounts[rowdex+curlib]+=temp[curlib]; }
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
			size_t n_to_drop=0;
			while (n_to_drop < ref_anchors.size() && ref_anchors[n_to_drop] < anchors[saved_dex] - bwidth) { 
                ++n_to_drop; 
            }
			if (n_to_drop) { 
				ref_anchors.erase(ref_anchors.begin(), ref_anchors.begin()+n_to_drop);
				ref_targets.erase(ref_targets.begin(), ref_targets.begin()+n_to_drop);
                for (int curlib=0; curlib<nlibs; ++curlib) { 
                    auto& cur_count=ref_count[curlib];
                    cur_count.erase(cur_count.begin(), cur_count.begin()+n_to_drop);
                }

				/* The 'matching_dex' gets pulled down as well. It's possible that the next 'saved_anchor'
				 * skips a couple of ref_anchors, so the number dropped might include 'matching_dex'. In
				 * that case, we just start from zero during the neighborhood calculations.
				 */
				ref_matching_dex-=n_to_drop;
				if (ref_matching_dex < 0) { 
                    ref_matching_dex=0; 
                }
			}
		} else if (engine.empty()) { break; }
	}

	// Storing into R objects.
    Rcpp::IntegerVector out_anchor1(anchors.begin(), anchors.end());
    for (auto& a1 : out_anchor1) { 
        a1+=fabin;
    }
    Rcpp::IntegerVector out_anchor2(targets.begin(), targets.end());
    for (auto& a2 : out_anchor2) { 
        a2 += ftbin;
    }

	const int ncombos=anchors.size();
    Rcpp::IntegerMatrix out_counts(ncombos, nlibs);
    auto cIt=counts.begin();
    for (int i=0; i<ncombos; ++i) {
        auto currow=out_counts.row(i);
        for (auto& val : currow) {
            val=*cIt;
            ++cIt;
        }
    }

    Rcpp::List out_ncounts(nmodes), out_areas(nmodes);
    for (size_t mode=0; mode<nmodes; ++mode) {
        const auto& ncounts=neighborcount[mode];
        const auto& narea=neighborarea[mode];
        Rcpp::IntegerMatrix curcounts(ncombos, nlibs);
        Rcpp::IntegerVector curarea(ncombos);
        
        if (mode >= startmode) { 
            auto ncIt=ncounts.begin();
            for (int i=0; i<ncombos; ++i) {
                auto ccrow=curcounts.row(i);
                for (auto& val : ccrow) {
                    val=*ncIt;
                    ++ncIt;
                }
            }
            std::copy(narea.begin(), narea.end(), curarea.begin());        
        }

        out_ncounts[mode]=curcounts;
        out_areas[mode]=curarea;
    }
      
	return Rcpp::List::create(out_anchor1, out_anchor2, out_counts, out_ncounts, out_areas);
    END_RCPP
}


