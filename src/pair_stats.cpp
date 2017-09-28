#include "diffhic.h"
#include "utils.h"

/* This function returns the fragment length, orientation and insert size 
 * corresponding to each processed read pair.
 */

SEXP pair_stats (SEXP anchor1_id, SEXP anchor2_id, SEXP anchor1_pos, SEXP anchor2_pos, SEXP anchor1_len, SEXP anchor2_len,
		SEXP same_chr, SEXP fstarts, SEXP fends) {
    BEGIN_RCPP

    const Rcpp::IntegerVector a1id(anchor1_id), a2id(anchor2_id), a1pos(anchor1_pos), 
        a2pos(anchor2_pos), a1len(anchor1_len), a2len(anchor2_len);
	const int np=LENGTH(a1id);
	if (np!=a2id.size() || np!=a1pos.size() || np!=a2pos.size() || np!=a1len.size() || np!=a2len.size()) {
		throw std::runtime_error("length of anchor/target position/length/index vectors must be equal"); 
	}

    const Rcpp::IntegerVector fs(fstarts), fe(fends);
	const int nf=fs.size();
	if (nf!=fe.size()) { throw std::runtime_error("length of fragment start and end vectors should be equal"); }

    const bool schr=check_logical_scalar(same_chr, "same chromosome specifier");
	
	// Running through the list of pairs.
    Rcpp::IntegerVector out_fraglen(np), out_ori(np), out_ins(np);
	
	for (int pair=0; pair<np; ++pair) {
	    int cural=a1len[pair];
		int curtl=a2len[pair];
		
		bool arev=cural < 0;
		bool trev=curtl < 0;
		if (arev) { cural *= -1; }
		if (trev) { curtl *= -1; }
		out_ori[pair] = (arev ? 1 : 0) + (trev ? 2 : 0);

		const int& curap=a1pos[pair];
		const int& curtp=a2pos[pair];
		const int curaend=curap+cural;
		const int curtend=curtp+curtl;
		if (!schr) {
			out_ins[pair]=NA_INTEGER;
		} else {
			out_ins[pair]=(curaend > curtend ? curaend : curtend) - (curap > curtp ? curtp : curap); 
			/* Compute insert size; protect against nested alignments, provide sensible results
			 * in cases where the apos of a reverse anchor read is below the tpos of a forward target read.
			 */
		}

 	    // Computing fragment lengths, unless fragment IDs are invalid.
		int& curflen=out_fraglen[pair];
        const int& aI=a1id[pair];
        const int& tI=a2id[pair];
        if (aI > 0 && tI > 0) { 
            if (aI > nf || tI > nf) {
                throw std::runtime_error("anchor indices out of range of fragments");
            }
            curflen += (arev ? curaend - fs[aI-1] : fe[aI-1] - curap + 1); // 1-based indexing.
            curflen += (trev ? curtend - fs[tI-1] : fe[tI-1] - curtp + 1);
        } else {
            curflen = NA_INTEGER;
        }
	}

	return Rcpp::List::create(out_fraglen, out_ori, out_ins);
    END_RCPP
}
