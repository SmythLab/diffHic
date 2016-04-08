#include "diffhic.h"
#include <cstdio>
#include <cstring>
#include "sam.h"

/***********************************************************************
 * Finds the fragment to which each read (or segment thereof) belongs.
 ***********************************************************************/

class base_finder {
public:
	base_finder() {}
    size_t nchrs() const { return pos.size(); }
	virtual int find_fragment(const int&, const int&, const bool&, const int&) const = 0;
    virtual ~base_finder() {};
protected:
	struct chr_stats {
		chr_stats(const int* s, const int* e, const int& l) : start_ptr(s), end_ptr(e), num(l) {}
		const int* start_ptr;
		const int* end_ptr;
		int num;
	};
	std::deque<chr_stats> pos;
};

class fragment_finder : public base_finder {
public:
	fragment_finder(SEXP, SEXP); // Takes a list of vectors of positions and reference names for those vectors.
	int find_fragment(const int&, const int&, const bool&, const int&) const;
};

fragment_finder::fragment_finder(SEXP starts, SEXP ends) {
	if (!isNewList(starts) || !isNewList(ends)) { throw std::runtime_error("start/end positions should each be a list of integer vectors"); }
    const int nchrs=LENGTH(starts);
	if (nchrs!=LENGTH(ends)) { throw std::runtime_error("number of start/end position vectors should be equal"); }
	
	for (int i=0; i<nchrs; ++i) {
		SEXP current1=VECTOR_ELT(starts, i);
		if (!isInteger(current1)) { throw std::runtime_error("start vector should be integer"); }
		SEXP current2=VECTOR_ELT(ends, i);
		if (!isInteger(current2)) { throw std::runtime_error("end vector should be integer"); }
		const int ncuts=LENGTH(current1);
		if (LENGTH(current2)!=ncuts) { throw std::runtime_error("start/end vectors should have the same length"); }
		pos.push_back(chr_stats(INTEGER(current1), INTEGER(current2), ncuts));	
	}
	return;
}

int fragment_finder::find_fragment(const int& c, const int& p, const bool& r, const int& l) const {
	// Binary search to obtain the fragment index with 5' end coordinates.
	int index=0;
	if (r) {
		int pos5=p+l-1;
		const int* eptr=pos[c].end_ptr;
		index= std::lower_bound(eptr, eptr+pos[c].num, pos5)-eptr;
		if (index==pos[c].num) {
			warning("read aligned off end of chromosome");
			--index;
		} else if (pos[c].start_ptr[index] > pos5) {
			warning("read aligned into spacer region for a filled-in genome");
			--index;
		}
	} else {
		/* If forward strand, we search relative to the start position of each fragment. If the
 		 * end of the identified fragment is less than the position, the read probably just
 		 * slipped below the start position of the true fragment and is partially sitting
 		 * in the spacer for filled-in genomes. We then simply kick up the index.
 		 */
		const int* sptr=pos[c].start_ptr;
		index=std::upper_bound(sptr, sptr+pos[c].num, p)-sptr-1;
		if (pos[c].end_ptr[index] < p) {
			warning("read starts in spacer region for a filled genome");
			++index;
		}
	}
	return index;
}

/***********************************************************************
 * Parses the CIGAR string to extract the alignment length, offset from 5' end of read.
 ***********************************************************************/

void parse_cigar (const bam1_t* read, int& alen, int& offset) {
    const uint32_t* cigar=bam_get_cigar(read);
    const int n_cigar=(read->core).n_cigar;
    if (n_cigar==0) {
        if ((read -> core).flag & BAM_FUNMAP) {
            alen=offset=0;
            return;
        }
        std::stringstream err;
        err << "zero-length CIGAR for read '" << bam_get_qname(read) << "'";
        throw std::runtime_error(err.str());
    }
	alen=bam_cigar2rlen(n_cigar, cigar);
    offset=0;

    if ((read->core).n_cigar) { 
        if (bam_is_rev(read)) { 
            const uint32_t last=cigar[(read->core).n_cigar - 1];
            if (bam_cigar_op(last)==BAM_CHARD_CLIP) {
                offset=bam_cigar_oplen(last);  
            }
        } else {
            const uint32_t first=cigar[0];
            if (bam_cigar_op(first)==BAM_CHARD_CLIP) {
                offset=bam_cigar_oplen(first);
            }
        }
    } 
	return;
}

/***********************************************
 * Something to hold segment, pair information.
 ***********************************************/

enum status { ISPET, ISMATE, NEITHER };

struct segment {
	int offset, alen, fragid, chrid, pos;
	bool reverse;
	bool operator<(const segment& rhs) const { return offset < rhs.offset; }
};

// These internal codes need to be negative, as fragment lengths are positive in return value.
#define internal_neither -1 
#define internal_ismate -2

int get_pet_dist (const segment& left, const segment& right) {
	if (right.chrid!=left.chrid || right.reverse==left.reverse) { return internal_neither; }
	const segment& fs=(left.reverse ? right : left);
	const segment& rs=(left.reverse ? left : right);
	if (fs.pos <= rs.pos) {
		if (fs.pos + fs.alen > rs.pos + rs.alen) { return internal_neither; }
		return rs.pos + rs.alen - fs.pos; // Assuming that alignment lengths are positive.
	}
	if (fs.pos < rs.pos+rs.alen) { return internal_neither; }
	return internal_ismate;
}

status get_status (const segment& left, const segment& right) {
	if (right.fragid!=left.fragid) { return NEITHER; }
	switch (get_pet_dist(left, right)) { 
		case internal_neither: return NEITHER;
		case internal_ismate: return ISMATE;
		default: return ISPET;
	}
}

/***********************************************
 * Something to identify invalidity of chimeric read pairs.
 ***********************************************/

struct check_invalid_chimera {
	virtual ~check_invalid_chimera() {};
	virtual bool operator()(const std::deque<segment>& read1, const std::deque<segment>& read2) const = 0;
};

struct check_invalid_by_fragid : public check_invalid_chimera {
	check_invalid_by_fragid() {};
	~check_invalid_by_fragid() {};
	bool operator()(const std::deque<segment>& read1, const std::deque<segment>& read2) const {
		if (read1.size()==2 && get_status(read2[0], read1[1])!=ISPET) { return true; }
		if (read2.size()==2 && get_status(read1[0], read2[1])!=ISPET) { return true; }
		return false;
	};
};

struct check_invalid_by_dist : public check_invalid_chimera {
	check_invalid_by_dist(SEXP span) {
		if (!isInteger(span) || LENGTH(span)!=1) { throw std::runtime_error("maximum chimeric span must be a positive integer"); }
		maxspan=asInteger(span);
	};
	~check_invalid_by_dist() {};
	bool operator()(const std::deque<segment>& read1, const std::deque<segment>& read2) const {
		if (read1.size()==2) {
			int temp=get_pet_dist(read2[0], read1[1]);
			if (temp < 0 || temp > maxspan) { return true; }
		}
		if (read2.size()==2) {
			int temp=get_pet_dist(read1[0], read2[1]);
			if (temp < 0 || temp > maxspan) { return true; }
		}
		/* Doesn't account for cases where the 5' end is nested inside the 3' end and the mate.
		 * These are physically impossible from a single linear DNA molecule. I suppose we can
		 * forgive this, because it could form from interactions betwen homologous chromosomes.
 		 */
		return false;
	};
	int get_span() const { return maxspan; }
private:
	int maxspan;
};

/************************
 * Something to read in BAM files; copied shamelessly from the 'bamsignals' package
 ************************/

class Bamfile {
public:
    Bamfile(const char * path) : holding(false) { 
        in = sam_open(path, "rb");
        if (in == NULL) { 
            std::stringstream out;
            out << "failed to open BAM file at '" << path << "'";
            throw std::runtime_error(out.str());
        }
        try {
            header=sam_hdr_read(in);
        } catch (std::exception &e) {
            sam_close(in);
            throw;
        }   
        read=bam_init1();
        next=bam_init1();
        return;
    }
    ~Bamfile(){
        sam_close(in);
        bam_hdr_destroy(header);
        bam_destroy1(read);
        bam_destroy1(next);
    }

    bool read_alignment() {
        if (holding) {
            bam1_t* tmp;            
            tmp=read;
            read=next;
            next=tmp;
            holding=false;
        } else {
            if (sam_read1(in, header, read) < 0) { return false; }
        }
        return true;
    }
       
    void put_back() {
        bam1_t* tmp;            
        tmp=read;
        read=next;
        next=tmp;
        holding=true;
        return;
    } 

    samFile* in;
    bam_hdr_t* header;
    bam1_t* read, *next;
    bool holding;
};

class OutputFile {
public: 
    OutputFile(const char* p, const int c1, const int c2) : out(NULL) { 
        std::stringstream converter;
        converter << p << c1 << "_" << c2;
        path=converter.str();
        return;
    }
    void add(int anchor, int target, int apos, int tpos, int alen, int tlen, bool arev, bool trev) {
        if (out==NULL) { 
            out = std::fopen(path.c_str(), "w");
            if (out==NULL) {
                std::stringstream err;
                err << "failed to open output file at '" << path << "'"; 
                throw std::runtime_error(err.str());
            }
        }
		if (alen<0 || tlen<0) { throw std::runtime_error("alignment lengths should be positive"); }
        if (arev) { alen *= -1; } 
        if (trev) { tlen *= -1; } 
        fprintf(out, "%i\t%i\t%i\t%i\t%i\t%i\n", anchor+1, target+1, apos, tpos, alen, tlen); // Get back to 1-indexing.
        return;
    }
    ~OutputFile() {
        if (out!=NULL) { std::fclose(out); }
    }
    
    FILE * out;
    std::string path;
};

/************************
 * Main loop.
 ************************/

SEXP internal_loop (const base_finder * const ffptr, status (*check_self_status)(const segment&, const segment&), const check_invalid_chimera * const icptr,
        SEXP chr_converter, SEXP bamfile, SEXP prefix, SEXP chimera_strict, SEXP minqual, SEXP do_dedup) {

    // Checking input values.
    if (!isString(bamfile) || LENGTH(bamfile)!=1) { throw std::runtime_error("BAM file path should be a character string"); }
    if (!isString(prefix) || LENGTH(prefix)!=1) { throw std::runtime_error("output file prefix should be a character string"); }
    if (!isLogical(chimera_strict) || LENGTH(chimera_strict)!=1) { throw std::runtime_error("chimera removal specification should be a logical scalar"); }
	if (!isLogical(do_dedup) || LENGTH(do_dedup)!=1) { throw std::runtime_error("duplicate removal specification should be a logical scalar"); }
	if (!isInteger(minqual) || LENGTH(minqual)!=1) { throw std::runtime_error("minimum mapping quality should be an integer scalar"); }

	// Initializing pointers.
    Bamfile input(CHAR(STRING_ELT(bamfile, 0)));
	const bool rm_invalid=asLogical(chimera_strict);
	const bool rm_dup=asLogical(do_dedup);
	const int minq=asInteger(minqual);
	const bool rm_min=!ISNA(minq);

    // Initializing the chromosome conversion table (to get from BAM TIDs to chromosome indices in the 'fragments' GRanges).
	const size_t nc=ffptr->nchrs();
    if (!isInteger(chr_converter)) { throw std::runtime_error("chromosome conversion table should be integer"); }
    if (LENGTH(chr_converter)!=nc) { throw std::runtime_error("conversion table should have length equal to the number of chromosomes"); }
    const int* converter=INTEGER(chr_converter);
    for (size_t i=0; i<nc; ++i) {
        if (converter[i]==NA_INTEGER || converter[i] < 0 || converter[i] >= nc) { throw std::runtime_error("conversion indices out of range"); }
    }
    
   	// Constructing output containers
    const char* oprefix=CHAR(STRING_ELT(prefix, 0));
	std::deque<std::deque<OutputFile> > collected(nc);
	for (size_t i=0; i<nc; ++i) { 
        for (size_t j=0; j<=i; ++j) { 
            collected[i].push_back(OutputFile(oprefix, i, j));
        }
    }
	int single=-1; // First one always reported as a singleton, as qname is empty.
	int total=0, dupped=0, filtered=0, mapped=0;
	int dangling=0, selfie=0;
	int total_chim=0, mapped_chim=0, multi_chim=0, inv_chimeras=0;

    bool isempty=false, isfirst=false, curdup=false, curunmap=false, hasfirst=false, hassecond=false,
         firstunmap=true, secondunmap=true, isdup=false, isunmap=false, ischimera=false;
    int nsegments=0;
   	std::deque<segment> read1, read2;
	segment current;
    std::string qname="";

    while (1) {
        isempty=true;
        nsegments=0;
        isfirst=curdup=curunmap=isdup=false;
        firstunmap=secondunmap=true;
        hasfirst=hassecond=false;
        read1.clear();
        read2.clear();

        while (input.read_alignment()) {
            isempty=false;
            if (std::strcmp(bam_get_qname(input.read), qname.c_str())!=0) { 
                // First one will pop out, but that's okay.
                qname=bam_get_qname(input.read);
                input.put_back();
                break;
            }
            ++nsegments;

            // Checking what the read is (first or second).
			isfirst=bool((input.read -> core).flag & BAM_FREAD1);
			if (isfirst) { hasfirst=true; }
			else { hassecond=true; }
            
			// Checking how we should proceed; whether we should bother adding it or not.
			curdup=bool((input.read -> core).flag & BAM_FDUP);
			curunmap=(bool((input.read -> core).flag & BAM_FUNMAP) || (rm_min && (input.read -> core).qual < minq));
			parse_cigar(input.read, current.alen, current.offset);
            if (current.offset==0 && current.alen > 0) { 
                if (curdup) { isdup=true; } // defaults to 'false' unless we have a definitive setting of markingness.
                if (!curunmap) { (isfirst ? firstunmap : secondunmap)=false; } // defaults to 'true' unless we know it's mapped (unmapped reads get alen=0 and won't reach here). 
            }

			// Checking which deque to put it in, if we're going to keep it.
            if (! (curdup && rm_dup) && ! curunmap) {
                current.reverse=bool(bam_is_rev(input.read));
                const int32_t& curtid=(input.read -> core).tid;
                if (curtid==-1 || curtid >= nc) {
                    std::stringstream err;
                    err << "tid for read '" << bam_get_qname(input.read) << "' out of range of BAM header";
                    throw std::runtime_error(err.str());
                } else {
                    current.chrid=converter[curtid];
                }
                current.pos=(input.read->core).pos + 1; // code assumes 1-based index for base position.

                std::deque<segment>& current_reads=(isfirst ? read1 : read2);
                if (current.offset==0) { current_reads.push_front(current); } 
                else { current_reads.push_back(current); }
            }
        }

        // If we processed nothing, then it's the end of the file and we break.
        if (isempty) { break; }

		// Skipping if it's a singleton; otherwise, reporting it as part of the total read pairs.
		if (!hasfirst || !hassecond) {
			++single;
			continue;
		}
		++total;

		// Adding to other statistics.
        ischimera=(nsegments > 2);
		if (ischimera) { ++total_chim; }
		if (isdup) { ++dupped; }
        isunmap=(firstunmap | secondunmap);
		if (isunmap) { ++filtered; }

		/* Skipping if unmapped, marked (and we're removing them), and if the first alignment
		 * of either read has any hard 5' clipping. This means that it's not truly 5' terminated
		 * (e.g. the actual 5' end was unmapped, duplicate removed or whatever). Note that
		 * not skipping UNMAP or DUP does not imply non-empty sets, as UNMAP/DUP are only set
		 * for 0-offset alignments; if this isn't in the file, these flags won't get set, but
		 * the sets can still be empty if non-zero-offset alignments are present and filtered
		 * (to escape the singles clause above). Thus, we need to check non-emptiness explicitly.
 		 */
		if (isunmap || (rm_dup && isdup) || read1.empty() || read2.empty() || read1.front().offset || read2.front().offset) { continue; }
		++mapped;

		// Assigning fragment IDs, if everything else is good.
		for (size_t i1=0; i1<read1.size(); ++i1) {
			segment& current=read1[i1];
			current.fragid=ffptr->find_fragment(current.chrid, current.pos, current.reverse, current.alen);
		}
		for (size_t i2=0; i2<read2.size(); ++i2) {
			segment& current=read2[i2];
			current.fragid=ffptr->find_fragment(current.chrid, current.pos, current.reverse, current.alen);
		}

		// Determining the type of construct if they have the same ID.
		switch ((*check_self_status)(read1.front(), read2.front())) {
			case ISPET:
				++dangling;
				continue;
			case ISMATE:
				++selfie;
				continue;
			default:
				break;
		}

		// Pulling out chimera diagnostics.
		if (ischimera) {
			++mapped_chim;
 		   	++multi_chim;	
			bool invalid=false;
			if (read1.size()==1 && read2.size()==1) {
				--multi_chim;
			} else if (read1.size() > 2 || read2.size() > 2) {
				invalid=true;
			} else {
				invalid=(*icptr)(read1, read2);
			}
			if (invalid) {
				++inv_chimeras;
				if (rm_invalid) { continue; }
			}
		}
		
		// Choosing the anchor segment, and reporting it.
		bool anchor=false;
		if (read1.front().chrid > read2.front().chrid) {
 		   anchor=true;
	   	} else if (read1.front().chrid==read2.front().chrid) {
			if (read1.front().fragid > read2.front().fragid) {
				anchor=true;
			} else if (read1.front().fragid == read2.front().fragid) {
				if (read1.front().pos > read2.front().pos) {
					anchor=true;
				}
			}
		}
		const segment& anchor_seg=(anchor ? read1.front() : read2.front());
		const segment& target_seg=(anchor ? read2.front() : read1.front());   
        collected[anchor_seg.chrid][target_seg.chrid].add(anchor_seg.fragid, target_seg.fragid, 
                anchor_seg.pos, target_seg.pos, anchor_seg.alen, target_seg.alen,
                anchor_seg.reverse, target_seg.reverse);
	}

	SEXP total_output=PROTECT(allocVector(VECSXP, 5));
	try {
        // Saving all file names.
        SET_VECTOR_ELT(total_output, 0, allocVector(VECSXP, nc));
        SEXP all_paths=VECTOR_ELT(total_output, 0);
        for (int i=0; i<nc; ++i) {
            SET_VECTOR_ELT(all_paths, i, allocVector(STRSXP, i+1));
            SEXP current_paths=VECTOR_ELT(all_paths, i);
            for (int j=0; j<=i; ++j) {
                if (collected[i][j].out!=NULL) {
                    SET_STRING_ELT(current_paths, j, mkChar(collected[i][j].path.c_str()));
                } else {
                    SET_STRING_ELT(current_paths, j, mkChar(""));
                }
            } 
        }

		// Dumping mapping diagnostics.
		SET_VECTOR_ELT(total_output, 1, allocVector(INTSXP, 4));
		int* dptr=INTEGER(VECTOR_ELT(total_output, 1));
		dptr[0]=total;
		dptr[1]=dupped;
		dptr[2]=filtered;
		dptr[3]=mapped;
	
		// Dumping the number of dangling ends, self-circles.	
		SET_VECTOR_ELT(total_output, 2, allocVector(INTSXP, 2));
		int * siptr=INTEGER(VECTOR_ELT(total_output, 2));
		siptr[0]=dangling;
		siptr[1]=selfie;

		// Dumping the number designated 'single', as there's no pairs.
		SET_VECTOR_ELT(total_output, 3, ScalarInteger(single));

		// Dumping chimeric diagnostics.
		SET_VECTOR_ELT(total_output, 4, allocVector(INTSXP, 4));
		int* cptr=INTEGER(VECTOR_ELT(total_output, 4));
		cptr[0]=total_chim;
		cptr[1]=mapped_chim;
		cptr[2]=multi_chim;
		cptr[3]=inv_chimeras;
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}
	UNPROTECT(1);
	return total_output;
}

SEXP report_hic_pairs (SEXP start_list, SEXP end_list, SEXP chrconvert, SEXP bamfile, SEXP outfile, 
        SEXP chimera_strict, SEXP chimera_span, SEXP minqual, SEXP do_dedup) try {
	fragment_finder ff(start_list, end_list);
	
	check_invalid_by_fragid invfrag; // Bit clunky to define both, but easiest to avoid nested try/catch.
	check_invalid_by_dist invdist(chimera_span);
	const check_invalid_chimera* invchim=NULL;
	if (invdist.get_span()==NA_INTEGER) { invchim=&invfrag; } 
	else { invchim=&invdist; }
	
	return internal_loop(&ff, &get_status, invchim, chrconvert, bamfile, outfile, chimera_strict, minqual, do_dedup);
} catch (std::exception& e) {
	return mkString(e.what());
}

/************************
 * Repeated loop for DNase Hi-C.
 ************************/

class simple_finder : public base_finder {
public:
	simple_finder(SEXP, SEXP);
	int find_fragment(const int&, const int&, const bool&, const int&) const;
private:
	int bin_width;
};

simple_finder::simple_finder(SEXP n_per_chr, SEXP bwidth) {
	if (!isInteger(bwidth)|| LENGTH(bwidth)!=1) { throw std::runtime_error("bin width must be an integer scalar"); }
	bin_width=asInteger(bwidth);

    if (!isInteger(n_per_chr)) { throw std::runtime_error("number of fragments per chromosome must be an integer vector"); }
    const int nchrs=LENGTH(n_per_chr);
    const int* nptr=INTEGER(n_per_chr);
	for (int i=0; i<nchrs; ++i) { pos.push_back(chr_stats(NULL, NULL, nptr[i])); }
	return;	
}

int simple_finder::find_fragment(const int& c, const int& p, const bool& r, const int& l) const {
	if (r) {
		int index=int((p+l-2)/bin_width); // -1, to get position of last base; -1 again, to get into the right bin.
		if (index >= pos[c].num) {
			warning("read aligned off end of chromosome");
			--index;
		}
		return index;
	}
	return int((p-1)/bin_width);
}

status no_status_check (const segment& left, const segment& right) {
	/* Fragment IDs have no concept in DNase Hi-C, so automatic
	 * detection of self-circles/dangling ends is impossible.
	 */
	(void)left;
	(void)right; // Just to avoid unused warnings, but maintain compatibility.
	return NEITHER;
}

SEXP report_hic_binned_pairs (SEXP num_in_chrs, SEXP bwidth, SEXP chrconvert, SEXP bamfile, SEXP outfile, 
        SEXP chimera_strict, SEXP chimera_span, SEXP minqual, SEXP do_dedup) try {
	simple_finder ff(num_in_chrs, bwidth);
	check_invalid_by_dist invchim(chimera_span);
	return internal_loop(&ff, &no_status_check, &invchim, chrconvert, bamfile, outfile, chimera_strict, minqual, do_dedup);
} catch (std::exception& e) {
	return mkString(e.what());
}

/********************
 * Testing functions.
 *******************/

SEXP test_parse_cigar (SEXP incoming) try {
	if (!isString(incoming) || LENGTH(incoming)!=1) { throw std::runtime_error("BAM file path should be a string"); }
    
    Bamfile input(CHAR(STRING_ELT(incoming, 0)));
    if (sam_read1(input.in, input.header, input.read)<0) { 
        throw std::runtime_error("BAM file is empty");
    } 

	SEXP output=PROTECT(allocVector(INTSXP, 2));
	int* optr=INTEGER(output);
	int& alen=*optr;
	int& offset=*(optr+1);
	parse_cigar(input.read, alen, offset);

    UNPROTECT(1);
	return(output);
} catch (std::exception& e) {
	return mkString(e.what());
}

SEXP test_fragment_assign(SEXP starts, SEXP ends, SEXP chrs, SEXP pos, SEXP rev, SEXP len) try {
	fragment_finder ff(starts, ends);
	if (!isInteger(chrs) || !isInteger(pos) || !isLogical(rev) || !isInteger(len)) { throw std::runtime_error("data types are wrong"); }
	const int n=LENGTH(chrs);
	if (n!=LENGTH(pos) || n!=LENGTH(rev) || n!=LENGTH(len)) { throw std::runtime_error("length of data vectors are not consistent"); }
	
	const int* cptr=INTEGER(chrs);
	const int* pptr=INTEGER(pos);
	const int* rptr=LOGICAL(rev);
	const int* lptr=INTEGER(len);

	SEXP output=PROTECT(allocVector(INTSXP, n));
	int *optr=INTEGER(output);
	for (int i=0; i<n; ++i) {
		optr[i]=ff.find_fragment(cptr[i], pptr[i], rptr[i], lptr[i])+1;
	}
	
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

