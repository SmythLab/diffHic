#include "diffhic.h"
#include "neighbors.h"
#include "utils.h"

SEXP quadrant_bg (SEXP anchor1, SEXP anchor2, SEXP count, 
		SEXP width, SEXP exclude, 
		SEXP anchor1_len, SEXP anchor2_len, SEXP issame) {
    BEGIN_RCPP

    const Rcpp::IntegerVector a1(anchor1), a2(anchor2), _count(count);
	const int npair=a1.size();
    if (npair!=a2.size() || npair!=_count.size()) { 
        throw std::runtime_error("input vectors must have the same length");
    }

    const int flank_width=check_integer_scalar(width, "flank width");
	const int exwidth=check_integer_scalar(exclude, "exclusion width");
    const int a1len=check_integer_scalar(anchor1_len, "first anchor length");
	const int a2len=check_integer_scalar(anchor2_len, "second anchor length");
    const bool intrachr=check_logical_scalar(issame, "same chromosome specifier");

    Rcpp::List out_counts(4), out_areas(4);
    for (int i=0; i<4; ++i) {
        out_counts[i]=Rcpp::IntegerVector(npair);
        out_areas[i]=Rcpp::IntegerVector(npair);
    }
	
	// Iterating over all quadrants.
	bottomright br(flank_width, a2len, intrachr, exwidth);		
	updown ud(flank_width, a2len, intrachr, exwidth);
	leftright lr(flank_width, a2len, intrachr, exwidth);
	allaround aa(flank_width, a2len, intrachr, exwidth);

	for (int quadtype=(intrachr ? 0 : 1); quadtype<4; ++quadtype) {
	    basic* current=NULL;
		switch(quadtype) { 
			case 0: current=&br; break;
			case 1: current=&ud; break;
			case 2: current=&lr; break;
			case 3: current=&aa; break;
		}

        Rcpp::IntegerVector ocounts=out_counts[quadtype]; // data passed by reference, so can still use for assignment.
        Rcpp::IntegerVector oarea=out_areas[quadtype];

		// Iterating across all flank widths.
		do {
			int running_sum=0;
			int left_index=0;
			int right_index=0;
			
			for (int curpair=0; curpair<npair; ++curpair) {
				current->set(a1[curpair], a2[curpair]);
				const int cur_anchor=current->row;
				if (cur_anchor >= a1len) { break; }
				const int left_edge=current->left;
				const int right_edge=current->right;

				// Identifying all bin pairs in the relevant range.
				while (left_index < npair && (a1[left_index] < cur_anchor || 
						(a1[left_index]==cur_anchor && a2[left_index] < left_edge))) {
					running_sum -= _count[left_index];
					++left_index;
				}

				while (right_index < npair && (a1[right_index]<cur_anchor || 
						(a1[right_index]==cur_anchor && a2[right_index] < right_edge))) { 
					running_sum += _count[right_index];
					++right_index;		
				}

                // Adding them onto the current location (except if the anchor is negative;
                // skipping needs to be done here, as running_sum still needs to be calculated).
                if (cur_anchor <  0) { continue; }
			    ocounts[curpair] += running_sum;
                oarea[curpair] += right_edge - left_edge;
			}
		} while (current->bump_level());
	}
	   
	return Rcpp::List::create(out_counts, out_areas);
    END_RCPP
}

