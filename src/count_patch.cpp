#include "read_count.h"
#include "utils.h"

SEXP count_patch(SEXP all, SEXP bin, SEXP filter, SEXP firstbin, SEXP lastbin) {
    BEGIN_RCPP
    const int f=check_integer_scalar(filter, "filter value");
    const int fbin=check_integer_scalar(firstbin, "index of first bin on second anchor chromosome");
	const int lbin=check_integer_scalar(lastbin, "index of last bin on second anchor chromosome");

	// Setting up the binning engine.
	binner engine(all, bin, fbin, lbin);
	const int nlibs=engine.get_nlibs();
	std::deque<int> counts, anchors, targets;

	// Running through all libraries.
	while (!engine.empty()) {
		engine.fill();
		const int curanchor=engine.get_anchor();
        const std::deque<int>& waschanged=engine.get_changed();
        const std::vector<int>& curcounts=engine.get_counts();

		// Adding it to the main list, if it's large enough.
        for (const auto& wc : waschanged) { 
			const int rowdex=wc*nlibs;
            const int enddex=rowdex+nlibs;
			const int countsum=std::accumulate(curcounts.begin()+rowdex, curcounts.begin()+enddex, 0);

			if (countsum >= f) { 
				anchors.push_back(curanchor);
				targets.push_back(wc + fbin);
				for (int index=rowdex; index<enddex; ++index) {
                    counts.push_back(curcounts[index]);
                }
			}
		}
	}

    Rcpp::IntegerVector out_a1(anchors.begin(), anchors.end());
    Rcpp::IntegerVector out_a2(targets.begin(), targets.end());
    const int ncombos=anchors.size();
    Rcpp::IntegerMatrix out_counts(ncombos, nlibs);
		
    // Iterating across and filling both the matrix and the components.
    auto cIt=counts.begin();
    for (size_t vecdex=0; vecdex<ncombos; ++vecdex) {
        auto currow=out_counts.row(vecdex);
        for (auto& val : currow) { 
            val=*cIt;
            ++cIt;
        }
    }

	return Rcpp::List::create(out_a1, out_a2, out_counts);
    END_RCPP
}
