#include "diffhic.h"

/* This function returns the fragment length, orientation and insert size 
 * corresponding to each processed read pair.
 */

SEXP pair_stats (SEXP aid, SEXP tid, SEXP apos, SEXP tpos, SEXP alen, SEXP tlen,
		SEXP same_chr, SEXP fstarts, SEXP fends) try {
	if (!isInteger(aid) || !isInteger(tid)) { throw std::runtime_error("anchor and target indices must be integer"); }
	if (!isInteger(apos) || !isInteger(tpos)) { throw std::runtime_error("anchor and target positions must be integer"); }
	if (!isInteger(alen) || !isInteger(tlen)) { throw std::runtime_error("anchor and target lengths must be integer"); }
	const int np=LENGTH(aid);
	if (np!=LENGTH(tid) || np!=LENGTH(apos) || np!=LENGTH(tpos) || np!=LENGTH(alen) || np!=LENGTH(tlen)) { 
		throw std::runtime_error("length of anchor/target position/length/index vectors must be equal"); 
	}

	if (!isLogical(same_chr) || LENGTH(same_chr)!=1) { throw std::runtime_error("same chromosome specifier should be a logical scalar"); }
	const bool schr=asLogical(same_chr);
	if (!isInteger(fstarts) || !isInteger(fends)) { throw std::runtime_error("fragment starts and ends should be integer vectors"); }
	const int nf=LENGTH(fstarts);
	if (nf!=LENGTH(fends)) { throw std::runtime_error("length of fragment start and end vectors should be equal"); }
	
	// Setting up pointers.
	const int* aiptr=INTEGER(aid),
		  *tiptr=INTEGER(tid),
		  *apptr=INTEGER(apos),
		  *tpptr=INTEGER(tpos),
		  *alptr=INTEGER(alen),
		  *tlptr=INTEGER(tlen);
	const int* fsptr=INTEGER(fstarts)-1,
		  *feptr=INTEGER(fends)-1;

	// Setting up output structures.
	SEXP output=PROTECT(allocVector(VECSXP, 3));
try {
	SET_VECTOR_ELT(output, 0, allocVector(INTSXP, np));
	SET_VECTOR_ELT(output, 1, allocVector(INTSXP, np));
	SET_VECTOR_ELT(output, 2, allocVector(INTSXP, np));
	int* foptr=INTEGER(VECTOR_ELT(output, 0)),
		* ooptr=INTEGER(VECTOR_ELT(output, 1)),
		* goptr=INTEGER(VECTOR_ELT(output, 2));

	// Setting up some objects.
	bool arev, trev;
	int cural, curtl, curaend, curtend;
	
	// Running through the list of pairs.
	for (int pair=0; pair<np; ++pair) {
		cural=alptr[pair];
		curtl=tlptr[pair];
		
		arev=cural < 0;
		trev=curtl < 0;
		if (arev) { cural *= -1; }
		if (trev) { curtl *= -1; }
		ooptr[pair] = (arev ? 1 : 0) + (trev ? 2 : 0);

		const int& curap=apptr[pair];
		const int& curtp=tpptr[pair];
		curaend=curap+cural;
		curtend=curtp+curtl;
		if (!schr) {
			goptr[pair]=NA_INTEGER;
		} else {
			goptr[pair]=(curaend > curtend ? curaend : curtend) - (curap > curtp ? curtp : curap); 
			/* Compute insert size; protect against nested alignments, provide sensible results
			 * in cases where the apos of a reverse anchor read is below the tpos of a forward target read.
			 */
		}

 	    // Computing fragment lengths.
		int& curflen=(foptr[pair]=0);
		curflen += (arev ? curaend - fsptr[aiptr[pair]] : feptr[aiptr[pair]] - curap + 1);
		curflen += (trev ? curtend - fsptr[tiptr[pair]] : feptr[tiptr[pair]] - curtp + 1);
	}
} catch (std::exception& e){
	UNPROTECT(1);
	throw;
}

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
