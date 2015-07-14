#include "diffhic.h"
#include "neighbours.h"

SEXP quadrant_bg (SEXP anchor, SEXP target, 
		SEXP abundance_int, SEXP abundance_dec, SEXP mult,
		SEXP width, SEXP exclude, 
		SEXP alen, SEXP tlen, SEXP issame) try {

	if (!isInteger(anchor) || !isInteger(target)) { throw std::runtime_error("anchor/target vectors must be integer"); }
	const int npair=LENGTH(anchor);
	if (LENGTH(target)!=npair) { throw std::runtime_error("anchor/target vectors must have the same length"); }
	if (!isInteger(abundance_int) || !isInteger(abundance_dec)) { throw std::runtime_error("vector of abundances should be integer"); }
	if (LENGTH(abundance_int)!=npair || LENGTH(abundance_dec)!=npair) { 
		throw std::runtime_error("vector of abundances should be the same length as that of the indices"); }
  	
	// Setting up pointers.
	const int * aptr=INTEGER(anchor), 
	  * tptr=INTEGER(target),
	  * biptr=INTEGER(abundance_int),
	  * bdptr=INTEGER(abundance_dec);

	// Determining the flank width.
	if (!isReal(mult) || LENGTH(mult)!=1) { throw std::runtime_error("multiplier must be a double-precision scalar"); }
	const double multiplier=asReal(mult);
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

	SEXP output=PROTECT(allocVector(REALSXP, npair));
	try {
		double * optr=REAL(output);
		int curpair;
		for (curpair=0; curpair<npair; ++curpair) { optr[curpair]=0; }
		int running_sum_int, running_sum_dec, 
			left_index, right_index,
			left_edge, right_edge, cur_anchor; 

		int* nptr=(int*)R_alloc(npair, sizeof(int));
		double* temp_int=(double*)R_alloc(npair, sizeof(double));
		double* temp_dec=(double*)R_alloc(npair, sizeof(double));
		double temp_val;

		// Iterating over all quadrants.
		bottomright br(flank_width, tlength, intrachr, exwidth);		
		updown ud(flank_width, tlength, intrachr, exwidth);
		leftright1 lr1(flank_width, tlength, intrachr, exwidth);
		leftright2 lr2(flank_width, tlength, intrachr, exwidth);
		allaround1 aa1(flank_width, tlength, intrachr, exwidth);
		allaround2 aa2(flank_width, tlength, intrachr, exwidth);
		basic* current=NULL;

		for (int quadtype=(intrachr ? 0 : 1); quadtype<6; ++quadtype) {
			switch(quadtype) { 
				case 0: current=&br; break;
				case 1: current=&ud; break;
				case 2: current=&lr1; break;
				case 3: current=&lr2; break;
				case 4: current=&aa1; break;
				case 5: current=&aa2; break;
			}

			if (quadtype!=3 && quadtype!=5) { // Don't clear them, we want them added on from 2 and 4, respectively.
				for (curpair=0; curpair<npair; ++curpair) { 
					nptr[curpair]=0; 
					temp_int[curpair]=0;
					temp_dec[curpair]=0;
				}
			}

			// Iterating across all flank widths.
			do {
				running_sum_int=0;
				running_sum_dec=0;
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
						running_sum_int -= biptr[left_index];
						running_sum_dec -= bdptr[left_index];
						++left_index;
					}

					while (right_index < npair && (aptr[right_index]<cur_anchor || 
							(aptr[right_index]==cur_anchor && tptr[right_index] < right_edge))) { 
						running_sum_int += biptr[right_index];
						running_sum_dec += bdptr[right_index];
						++right_index;		
					}

					if (cur_anchor >= 0) {
						temp_int[curpair] += running_sum_int;
						temp_dec[curpair] += running_sum_dec;
						nptr[curpair] += right_edge - left_edge; // Figuring out the actual number of boxes.
					}
				}
			} while (current->bump_level());

			// Checking if it exceeds the previous maxima (2 and 4 are only half done, so we skip them).
			if (quadtype!=2 && quadtype!=4) { 
				for (curpair=0; curpair<npair; ++curpair) {
					if (nptr[curpair]) {
						temp_val = (temp_int[curpair] + temp_dec[curpair]/multiplier)/nptr[curpair];
						if (optr[curpair] < temp_val) { optr[curpair]=temp_val; }
//						if (exwidth) { Rprintf("%i %i %.3f %i\n", aptr[curpair]+1, tptr[curpair]+1, temp_val, nptr[curpair]); }
					}
				}
			}
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


