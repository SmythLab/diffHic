#include "read_count.h"
#include "utils.h"

SEXP directionality(SEXP all, SEXP bin, SEXP span, SEXP first_bin, SEXP last_bin) {
    BEGIN_RCPP

	// Getting scalar values.
    const int fbin=check_integer_scalar(first_bin, "index of first bin");
    const int lbin=check_integer_scalar(last_bin, "index of last bin");
    const size_t sp=check_integer_scalar(span, "span to compute directionality");

	// Setting up the binning engine.
	binner engine(all, bin, fbin, lbin);
	const int nlibs=engine.get_nlibs();
	const int nbins=engine.get_nbins();

    // Settin gup the memory containers.
	std::deque<int> counts, anchors, targets;

	// Setting up the output vectors immediately.
    Rcpp::IntegerMatrix outdown(nbins, nlibs);
    Rcpp::IntegerMatrix outup(nbins, nlibs);

	while (!engine.empty()) {
		engine.fill();
		int curanchor=engine.get_anchor() - fbin;
        const std::deque<int>& waschanged=engine.get_changed();
        const std::vector<int>& curcounts=engine.get_counts();

        for (const auto& rowdex : waschanged) { 
			size_t diff=curanchor-rowdex;

			// Filling up the directionality indices.		
			if (diff && diff <= sp) { 
                auto downrow=outdown.row(curanchor);
                auto dIt=downrow.begin();
                auto uprow=outup.row(rowdex);
                auto uIt=uprow.begin();

                auto ccIt=curcounts.begin() + rowdex*nlibs;
                for (int lib=0; lib<nlibs; ++lib, ++ccIt, ++dIt, ++uIt) {
                    const int& thiscount=(*ccIt);
                    (*dIt)+=thiscount;
                    (*uIt)+=thiscount;
                }
			} 
		}
	}

	return Rcpp::List::create(outdown, outup);
    END_RCPP
}

