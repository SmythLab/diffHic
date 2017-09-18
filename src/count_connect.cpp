#include "diffhic.h"
#include "read_count.h"
#include "utils.h"

void accumulate_pairs(int curab, int curtb, pair_queue& next, 
        const std::vector<Rcpp::IntegerVector>& anchor1, const std::vector<Rcpp::IntegerVector>& anchor2,
        const std::vector<int>& nums, std::vector<int>& indices, std::vector<int>& counts) {
    
    // This speeds things up by collecting entries with the same fragment indices and determining its count vector.
    std::fill(counts.begin(), counts.end(), 0);
    do {
        int curlib=next.top().library;
        int& libdex=indices[curlib];
        counts[curlib]+=1;
        next.pop();
        if ((++libdex) < nums[curlib]) {
            next.push(coord(anchor1[curlib][libdex], anchor2[curlib][libdex], curlib)); // 1-based indexing.
        } 
    } while (!next.empty() && next.top().anchor==curab && next.top().target==curtb);
    return;
}

SEXP count_connect(SEXP all, SEXP start1, SEXP end1, SEXP region1, 
        SEXP start2, SEXP end2, SEXP region2, SEXP filter) {		
    BEGIN_RCPP

    Rcpp::IntegerVector _start1(start1), _end1(end1);
    const int ni1=_start1.size();
    if (_end1.size()!=ni1) { 
        throw std::runtime_error("start/end index vectors (1) should be the same length"); 
    }
    auto s1It=_start1.begin()-1, e1It=_end1.begin()-1; // Account for 1-based indexing.

    Rcpp::IntegerVector _region1(region1);
    const int nrp1=_region1.size()+1; // As 'end' refers to one-past-the-end. 
    auto r1It=_region1.begin()-1;
    
    const int filtval=check_integer_scalar(filter, "filter value");

    // Repeating for the second region, if specified.
    Rcpp::IntegerVector _start2(_start1), _end2(_end1), _region2(_region1);
    int ni2=ni1, nrp2=nrp1;

    const bool use_second=(region2!=R_NilValue);
    if (use_second) {
        _start2=Rcpp::IntegerVector(start2);
        _end2=Rcpp::IntegerVector(end2);
        ni2=_start2.size();
        if (_end2.size()!=ni2) { 
            throw std::runtime_error("start/end index vectors (2) should be the same length"); 
        }

        _region2=Rcpp::IntegerVector(region2);
        nrp2=_region2.size()+1;
    }
    
    auto s2It=_start2.begin()-1, e2It=_end2.begin()-1, r2It=_region2.begin()-1;

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
	std::vector<Rcpp::IntegerVector> anchor1, anchor2;
	std::vector<int> nums, indices;
    const int nlibs=setup_pair_data(all, anchor1, anchor2, nums, indices);
    
    pair_queue next;
    for (int lib=0; lib<nlibs; ++lib) {
        if (nums[lib]) { next.push(coord(anchor1[lib][0], anchor2[lib][0], lib)); } // 1-based indexing.
    }
	
	// Running through all libraries.
	std::deque<int> counts;
	typedef std::pair<int, int> combo;
	std::map<combo, std::pair<int, int> > bins;
	std::vector<int> curcounts(nlibs);
    int counter=0;

	int smallest=-1;
	std::deque<int> returned_anchors, returned_targets;
	std::deque<size_t> count_indices;

	while (!next.empty()) {
		const int curab=next.top().anchor;
		const int curtb=next.top().target;
		++counter;
        accumulate_pairs(curab, curtb, next, anchor1, anchor2, nums, indices, curcounts); 

		/* Allocating counts to every pair of ranges containing these fragments. This
         * has a lot of loops but it shouldn't be too bad if there aren't many 
         * overlapping ranges. The 'uppermode' just checks the anchor vs the second
         * set of regions and the target vs the first set of regions, if there are two sets.
 		 */
        for (int mode=0; mode<(use_second ? 2 : 1); ++mode) { 
	        int y1, y2;
            if (mode==0) { 
                y1=curab;
                y2=curtb;
            } else {
                y1=curtb;
                y2=curab;
            }
            if (y1>ni1) { throw std::runtime_error("invalid index (1) for supplied fragments"); } // 1-based indexing, so '>' is right.
            const int& s1x=*(s1It + y1);
            const int& e1x=*(e1It + y1);
            if (y2>ni2) { throw std::runtime_error("invalid index (2) for supplied fragments"); }
            const int& s2x=*(s2It + y2);
            const int& e2x=*(e2It + y2);
        
            if (s1x!=e1x && s2x!=e2x) { 
                if (s1x <= 0 || s2x <= 0 || e1x > nrp1 || e2x > nrp2) { throw std::runtime_error("invalid start/endpoints for region indices"); }
                for (int x1=s1x; x1<e1x; ++x1) {
                    for (int x2=s2x; x2<e2x; ++x2) {

                        /* Avoid redundant naming, when looking within the same ranges; region with 
                         * higher genomic coordinates goes first. This should always increase, see below.
                         */
	                    combo temp;
                        const int& r1=*(r1It + x1);
                        const int& r2=*(r2It + x2);
                        if (r1 > r2) { 
                            temp.first=r1;
                            temp.second=r2;
                        } else {
                            temp.first=r2;
                            temp.second=r1;
    					}
	    				if (temp.first < smallest) { throw std::runtime_error("largest index cannot be lower than 'smallest' threshold"); }

    					auto itb=bins.lower_bound(temp);
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
    					for (int lib=0; lib<nlibs; ++lib) { counts[index+lib] += curcounts[lib]; }
                    }
				}
			}
		}

		/* Saving all entries where the first region index is below the smallest anchor-overlapped region index for the current read pair.
         * This frees up the map to make it almost 1-dimensional, which should reduce memory usage and improve speed.
         *
         * The logic is that, as the anchor index increases, the anchor fragment cannot overlap a region with a lower index in 'r(1/2)ptr'
         * (assuming that the regions have been sorted at the R level). So, the anchor-overlapped region index cannot decrease.
         * Now, the first index of each entry is the larger of the two indices; if an entry has a first index below the current
         * smallest anchor-overlapped region, it cannot possibly be updated at later iterations. This means we can remove it.
  		 */
        const int& s1x=*(s1It+curab);
        const int& e1x=*(e1It+curab);
        const int altsmallest1=*std::min_element(r1It+s1x, r1It+e1x);
        if (altsmallest1 < smallest) { 
            smallest=altsmallest1;
        }
        if (use_second) { 
            const int& s2x=*(s2It+curab);
            const int& e2x=*(e2It+curab);
            const int altsmallest2=*std::min_element(r2It+s2x, r2It+e2x);
            if (altsmallest2 < smallest) { 
                smallest=altsmallest2;
            }
        }

		auto itb=bins.begin();
		while (itb!=bins.end() && (itb->first).first < smallest) {
            const int& index=(itb->second).first;
            const int countsum=std::accumulate(counts.begin()+index, counts.begin()+index+nlibs, 0);
    
            if (countsum >= filtval) { // Checking that the count sum is sufficient.
                returned_anchors.push_back((itb->first).first);
                returned_targets.push_back((itb->first).second);
                count_indices.push_back((itb->second).first);
            }
            bins.erase(itb++);
		}
	}

	// Assessing how many remaining combinations are above threshold.
	for (auto itb=bins.begin(); itb!=bins.end(); ++itb) { 
        const int& index=(itb->second).first;
        const int countsum=std::accumulate(counts.begin()+index, counts.begin()+index+nlibs, 0);
		if (countsum >= filtval) { 
			returned_anchors.push_back((itb->first).first);
			returned_targets.push_back((itb->first).second);
			count_indices.push_back((itb->second).first);
		}
	}

	// Returning all count combinations underneath the threshold.
    const size_t ncombos=count_indices.size();
    Rcpp::IntegerVector out_anchor1(returned_anchors.begin(), returned_anchors.end());
    Rcpp::IntegerVector out_anchor2(returned_targets.begin(), returned_targets.end());
    Rcpp::IntegerMatrix out_count(ncombos, nlibs);
		
	for (size_t odex=0; odex < ncombos; ++odex) { 
		int index=count_indices[odex];
        auto currow=out_count.row(odex);
		for (auto& val : currow) {
            val=counts[index];
            ++index;
        }
    }

	return Rcpp::List::create(out_anchor1, out_anchor2, out_count);
    END_RCPP
}

/* This does the same for DNase-C data, where the assignment into ranges 
 * has already been done at the R-level using linkOverlaps. So we just
 * need to scan across the assignments in 'links' and aggregate counts.
 */

SEXP count_reconnect(SEXP links, SEXP filter) {	
    BEGIN_RCPP 

	// Setting up the structure to the links.
	std::vector<Rcpp::IntegerVector> anchor1, anchor2;
	std::vector<int> nums, indices;
 	const int nlibs=setup_pair_data(links, anchor1, anchor2, nums, indices);

    pair_queue next;
    for (int lib=0; lib<nlibs; ++lib) {
        if (nums[lib]) { next.push(coord(anchor1[lib][0], anchor2[lib][0], lib)); }
    }

    const int filtval=check_integer_scalar(filter, "filter value");

   // Running through the pairs.
    std::deque<int> returned_anchors, returned_targets;
	std::deque<int> counts;
	std::vector<int> curcounts(nlibs);

    while (!next.empty()) {
		const int curab=next.top().anchor;
		const int curtb=next.top().target;
        accumulate_pairs(curab, curtb, next, anchor1, anchor2, nums, indices, curcounts); 

        // Inserting if the sum of counts exceeds the value.
        const int total=std::accumulate(curcounts.begin(), curcounts.end(), 0);
        if (total >= filtval) { 
            returned_anchors.push_back(curab);
            returned_targets.push_back(curtb);
            counts.insert(counts.end(), curcounts.begin(), curcounts.end());
        }
    }

    // Saving the output.
    Rcpp::IntegerVector out_a1(returned_anchors.begin(), returned_anchors.end());
    Rcpp::IntegerVector out_a2(returned_targets.begin(), returned_targets.end());
    const int ncombos = returned_anchors.size();
    Rcpp::IntegerMatrix out_counts(ncombos, nlibs);

    auto cIt=counts.begin();
    for (size_t odex=0; odex < ncombos; ++odex) { 
        auto currow=out_counts.row(odex);
		for (auto& val : currow) { 
            val=*cIt;
            ++cIt;
        }
	}

    return Rcpp::List::create(out_a1, out_a2, out_counts);
    END_RCPP
}

