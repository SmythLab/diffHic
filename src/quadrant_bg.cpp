#include "diffhic.h"
#include "neighbors.h"

SEXP quadrant_bg (SEXP anchor, SEXP target, SEXP count, 
		SEXP width, SEXP exclude, 
		SEXP alen, SEXP tlen, SEXP issame) try {

	if (!isInteger(anchor) || !isInteger(target)) { throw std::runtime_error("anchor/target vectors must be integer"); }
	const int npair=LENGTH(anchor);
	if (LENGTH(target)!=npair) { throw std::runtime_error("anchor/target vectors must have the same length"); }
	if (!isInteger(count)) { throw std::runtime_error("vector of abundances should be integer"); }
	if (LENGTH(count)!=npair) { throw std::runtime_error("vector of counts should be the same length as that of the indices"); }
  	
	// Setting up pointers.
	const int * aptr=INTEGER(anchor), 
	  * tptr=INTEGER(target),
      * cptr=INTEGER(count);

	// Determining the flank width.
	if (!isInteger(width) || LENGTH(width)!=1) { throw std::runtime_error("flank width must be an integer scalar"); }
	const int flank_width=asInteger(width);
	if (!isInteger(exclude) || LENGTH(exclude)!=1) { throw std::runtime_error("exclusion width must be an integer scalar"); }
	const int exwidth=asInteger(exclude);

	if (!isInteger(alen) || LENGTH(alen)!=1) { throw std::runtime_error("anchor length must be an integer scalar"); }
	const int alength=asInteger(alen);
	if (!isInteger(tlen) || LENGTH(tlen)!=1) { throw std::runtime_error("anchor length must be an integer scalar"); }
	const int tlength=asInteger(tlen);
	if (!isLogical(issame) || LENGTH(issame)!=1) { throw std::runtime_error("same chromosome specifier must be a logical scalar"); }
	const bool intrachr=asLogical(issame);

	SEXP output=PROTECT(allocVector(VECSXP, 2));
	try {
        // Setting up output for the counts.
        SET_VECTOR_ELT(output, 0, allocVector(VECSXP, 4));
        SEXP count_out=VECTOR_ELT(output, 0); 
        int** optrs=(int**)R_alloc(4, sizeof(int*));
        for (int i=0; i<4; ++i) {
            SET_VECTOR_ELT(count_out, i, allocVector(INTSXP, npair));
            optrs[i]=INTEGER(VECTOR_ELT(count_out, i));
            std::fill(optrs[i], optrs[i]+npair, 0);
        }

        // Setting up output for the neighbourhood size.
        SET_VECTOR_ELT(output, 1, allocVector(VECSXP, 4));
        SEXP n_out=VECTOR_ELT(output, 1); 
        int** nptrs=(int**)R_alloc(4, sizeof(int*));
        for (int i=0; i<4; ++i) {
            SET_VECTOR_ELT(n_out, i, allocVector(INTSXP, npair));
            nptrs[i]=INTEGER(VECTOR_ELT(n_out, i));
            std::fill(nptrs[i], nptrs[i]+npair, 0);
        }

        int* optr=NULL, *nptr=NULL;
		int curpair, running_sum,
			left_index, right_index,
			left_edge, right_edge, cur_anchor; 

		// Iterating over all quadrants.
		bottomright br(flank_width, tlength, intrachr, exwidth);		
		updown ud(flank_width, tlength, intrachr, exwidth);
		leftright lr(flank_width, tlength, intrachr, exwidth);
		allaround aa(flank_width, tlength, intrachr, exwidth);
		basic* current=NULL;

		for (int quadtype=(intrachr ? 0 : 1); quadtype<4; ++quadtype) {
			switch(quadtype) { 
				case 0: current=&br; break;
				case 1: current=&ud; break;
				case 2: current=&lr; break;
				case 3: current=&aa; break;
			}
            optr=optrs[quadtype];
            nptr=nptrs[quadtype];

			// Iterating across all flank widths.
			do {
				running_sum=0;
				left_index=0;
				right_index=0;
				
				for (curpair=0; curpair<npair; ++curpair) {
					current->set(aptr[curpair], tptr[curpair]);
					cur_anchor=current->row;
					if (cur_anchor >= alength) { break; }
					left_edge=current->left;
					right_edge=current->right;

					// Identifying all bin pairs in the relevant range.
					while (left_index < npair && (aptr[left_index] < cur_anchor || 
							(aptr[left_index]==cur_anchor && tptr[left_index] < left_edge))) {
						running_sum -= cptr[left_index];
						++left_index;
					}

					while (right_index < npair && (aptr[right_index]<cur_anchor || 
							(aptr[right_index]==cur_anchor && tptr[right_index] < right_edge))) { 
						running_sum += cptr[right_index];
						++right_index;		
					}

                    // Adding them onto the current location (except if the anchor is negative;
                    // skipping needs to be done here, as running_sum still needs to be calculated).
                    if (cur_anchor <  0) { continue; }
				    optr[curpair] += running_sum;
                    nptr[curpair] += right_edge - left_edge;
				}
			} while (current->bump_level());
		}
	} catch (std::exception &e) {
		UNPROTECT(1);
		throw;
	}
	   
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}


