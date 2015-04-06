#include "diffhic.h"

/* Setting up a couple of structures to direct how the additions are performed.
 */

struct basic {
	basic(int w, int t, bool i) : level(0), width(w), tlen(t), intra(i), remove_self(false) {}
	virtual void set (int, int)=0;
	virtual ~basic() {};
	virtual bool bump_level () { 
		if (level >= width) { return false; }
		++level;
		return true;
	}
	bool discard_self() { return (level==0 && remove_self); }
	int row, left, right;
protected:
	int level, width, tlen;
	bool intra, remove_self;

	void restrain () {
		if (left < 0) { left=0; }
		if (intra) {
			if (right > row) { 
				right=row+1; 
				if (left > right) { left=right; }
			}
		} else if (right > tlen) { right=tlen; } // For intra's, right will hit diagonal; no need to worry about tlen.
	}
};

struct bottomright : public basic { 
	bottomright(int w, int t, bool i) : basic(w, t, i) {}
	~bottomright() {};
	void set(int a, int t) {
		row=a-level;
		left=(level ? t : t+1); 
		right=t+width+1; 
		restrain();
	}
};

struct updown : public basic {
	updown(int w, int t, bool i) : basic(w, t, i) { level=-w; }
	~updown() {};
	void set(int a, int t) {
		row=a+level;
		left=t;
		right=(level==0 ? t : t+1); 
		restrain();
	}	
};

struct leftright : public basic {
	leftright(int w, int t, bool i) : basic(w, t, i) { remove_self=true; }
	~leftright() {};
	bool bump_level() { return false; }
	void set(int a, int t) { 
		row=a;
		left=t-width;
		right=t+width+1; 
		restrain(); 
	}	
};

struct allaround : public basic {
	allaround(int w, int t, bool i) : basic(w, t, i) { 
		remove_self=true; 
		level=-w;
	}
	~allaround() {};
	void set(int a, int t) {
		row=a+level;
		left=t-width;
		right=t+width+1;
		restrain();		
	}
};

/* Main loop */

extern "C" {

SEXP quadrant_bg (SEXP anchor, SEXP target, 
		SEXP abundance_int, SEXP abundance_dec, SEXP mult,
		SEXP width, SEXP alen, SEXP tlen, SEXP issame) try {
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
		bottomright br(flank_width, tlength, intrachr);		
		updown ud(flank_width, tlength, intrachr);
		leftright lr(flank_width, tlength, intrachr);
		allaround aa(flank_width, tlength, intrachr);
		basic* current=NULL;

		for (int quadtype=(intrachr ? 0 : 1); quadtype<4; ++quadtype) {
			switch(quadtype) { 
				case 0: current=&br; break;
				case 1: current=&ud; break;
				case 2: current=&lr; break;
				case 3: current=&aa; break;
			}
			for (curpair=0; curpair<npair; ++curpair) { 
				nptr[curpair]=0; 
				temp_int[curpair]=0;
				temp_dec[curpair]=0;
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
						if (!current->discard_self()) { 
							temp_int[curpair] += running_sum_int;
							temp_dec[curpair] += running_sum_dec;
							nptr[curpair] += right_edge - left_edge; // Figuring out the actual number of boxes.
						} else {
							temp_int[curpair] += running_sum_int - biptr[curpair]; 
							temp_dec[curpair] += running_sum_dec - bdptr[curpair];
							nptr[curpair] += right_edge - left_edge - 1;
//						Rprintf("%i L:%i R:%i %i %i\n", cur_anchor, left_edge, right_edge, aptr[curpair], tptr[curpair]);
						}

					}
				}
			} while (current->bump_level());

			// Checking if it exceeds the previous maxima.
			for (curpair=0; curpair<npair; ++curpair) {
				if (nptr[curpair]) {
					temp_val = (temp_int[curpair] + temp_dec[curpair]/multiplier)/nptr[curpair];
					if (optr[curpair] < temp_val) { optr[curpair]=temp_val; }
//	 			    Rprintf("%i %i %.3f\n", aptr[curpair]+1, tptr[curpair]+1, temp_val);
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

}
