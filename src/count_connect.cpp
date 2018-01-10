#include "diffhic.h"
#include "read_count.h"
#include "utils.h"

SEXP count_connect(SEXP links, SEXP filter) {	
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
    
        // This speeds things up by collecting entries with the same fragment indices and determining its count vector.
        std::fill(curcounts.begin(), curcounts.end(), 0);
        do {
            int curlib=next.top().library;
            int& libdex=indices[curlib];
            curcounts[curlib]+=1;
            next.pop();
            if ((++libdex) < nums[curlib]) {
                next.push(coord(anchor1[curlib][libdex], anchor2[curlib][libdex], curlib)); // 1-based indexing.
            } 
        } while (!next.empty() && next.top().anchor==curab && next.top().target==curtb);

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

