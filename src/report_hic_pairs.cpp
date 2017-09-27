#include "diffhic.h"
#include "utils.h"

#include <cstdio>
#include <cstring>
#include "sam.h"

/***********************************************
 * Something to hold segment, pair information.
 ***********************************************/

struct segment {
    segment() : offset(0), width(0), chrid(0), pos(0), fragid(NA_INTEGER), reverse(false) {}
    segment(int c, int p, bool r, int o, int w) : offset(o), width(w), chrid(c), pos(p), fragid(NA_INTEGER), reverse(r) {}
	const int offset, width, chrid, pos;
    int fragid;
	const bool reverse;
    int get_5pos() const { return (reverse ? pos + width - 1 : pos); }
};

/***********************************************************************
 * Finds the fragment to which each read (or segment thereof) belongs.
 ***********************************************************************/

class base_finder {
public:
	base_finder() {}
    virtual size_t nchrs() const=0;
	virtual int find_fragment(const segment&) const=0;
    virtual ~base_finder() {};
};

// A class for typical Hi-C experiments.

class fragment_finder : public base_finder {
public:
    fragment_finder(SEXP, SEXP);
    size_t nchrs() const;
	int find_fragment(const segment&) const;
private:
	std::vector<Rcpp::IntegerVector> fstarts, fends;
};

size_t fragment_finder::nchrs() const {
    return fstarts.size();
}

fragment_finder::fragment_finder(SEXP starts, SEXP ends) { // Takes a list of vectors of start/end fragment positions for each chromosome.
    Rcpp::List _starts(starts), _ends(ends);
    const int nchrs=_starts.size();
	if (nchrs!=_ends.size()) { throw std::runtime_error("number of start/end position vectors should be equal"); }
    fstarts.resize(nchrs);
    fends.resize(nchrs);

	for (int i=0; i<nchrs; ++i) {
        Rcpp::IntegerVector curstarts=_starts[i], curends=_ends[i];
		const int ncuts=curstarts.size();
		if (curends.size()!=ncuts) { 
            throw std::runtime_error("start/end vectors should have the same length"); 
        }
        fstarts[i]=curstarts;
        fends[i]=curends;
	}
	return;
}

int fragment_finder::find_fragment(const segment& current) const {
    const int& c=current.chrid;
    const int pos5=current.get_5pos();
        
	// Binary search to obtain the fragment index with 5' end coordinates.
	if (current.reverse) {
		const auto& fend=fends[c];
        int index=std::lower_bound(fend.begin(), fend.end(), pos5)-fend.begin();
        if (index==fend.size()) { 
            Rcpp::warning("read aligned off end of chromosome");
			--index;
		}
        return index;
	} else {
        const auto& fstart=fstarts[c];
		return (std::upper_bound(fstart.begin(), fstart.end(), pos5)-fstart.begin())-1;
	}
}

/***********************************************************************
 * Parses the CIGAR string to extract the alignment length, offset from 5' end of read.
 ***********************************************************************/

void parse_cigar (const bam1_t* read, int& offset, int& width) {
    const uint32_t* cigar=bam_get_cigar(read);
    const int n_cigar=(read->core).n_cigar;
    if (n_cigar==0) {
        if ((read -> core).flag & BAM_FUNMAP) {
            width=offset=0;
            return;
        }
        std::stringstream err;
        err << "zero-length CIGAR for read '" << bam_get_qname(read) << "'";
        throw std::runtime_error(err.str());
    }
	width=bam_cigar2rlen(n_cigar, cigar);
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

/***********************************************************************
 * Determine if two segments are paired-end and inward-facing.
 ***********************************************************************/

enum status { ISPET, ISMATE, NEITHER };

int get_pet_dist (const segment& left, const segment& right, status& flag) { 
    /* Computing distances between 5' ends. Overextended and nested reads are just truncated here,
     * attributable to trimming failures, as they are impossible to generate from a single DNA
     * fragment or ligation product. We also assume that alignment lengths are positive. 
     */
	if (right.chrid!=left.chrid || right.reverse==left.reverse) { 
        flag=NEITHER;
        return 0;
    }
    int f5, r5;
    if (left.reverse) {
        f5=right.pos;
        r5=left.get_5pos();
    } else {
        f5=left.pos;
        r5=right.get_5pos();
    }
    if (r5 < f5) { 
        flag=ISMATE;
        return 0;
    }
    flag=ISPET;
    return r5 - f5 + 1;
}

status get_status (const segment& left, const segment& right) {
	if (right.fragid!=left.fragid) { return NEITHER; }
    status flag;
    get_pet_dist(left, right, flag);
    return flag;
}

/***********************************************
 * Something to identify invalidity of chimeric read pairs.
 ***********************************************/

struct check_invalid_chimera { // virtual class
	virtual ~check_invalid_chimera() {};
	virtual bool operator()(const std::deque<segment>& read1, const std::deque<segment>& read2) const = 0;
};

struct check_invalid_by_fragid : public check_invalid_chimera { // check based on fragment ID.
	check_invalid_by_fragid() {};
	~check_invalid_by_fragid() {};
	bool operator()(const std::deque<segment>& read1, const std::deque<segment>& read2) const {
		if (read1.size()==2 && get_status(read2[0], read1[1])!=ISPET) { return true; }
		if (read2.size()==2 && get_status(read1[0], read2[1])!=ISPET) { return true; }
		return false;
	};
};

struct check_invalid_by_dist : public check_invalid_chimera { // check based on distance.
	check_invalid_by_dist(SEXP span) : maxspan(check_integer_scalar(span, "maximum chimeric span")) {}

	~check_invalid_by_dist() {};

	bool operator()(const std::deque<segment>& read1, const std::deque<segment>& read2) const {
        status flag;
        int temp;
		if (read1.size()==2) {
			temp=get_pet_dist(read2[0], read1[1], flag);
			if (flag!=ISPET || temp > maxspan) { return true; }
		}
		if (read2.size()==2) {
			temp=get_pet_dist(read1[0], read2[1], flag);
			if (flag!=ISPET || temp > maxspan) { return true; }
		}
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
    OutputFile(const char* p, const int c1, const int c2, const size_t np) : num(0), NPAIRS(np), 
            ai(NPAIRS), ti(NPAIRS), ap(NPAIRS), tp(NPAIRS), al(NPAIRS), tl(NPAIRS), out(NULL), saved(false) {
        std::stringstream converter;
        converter << p << c1 << "_" << c2;
        path=converter.str();
    }

    void add(const segment& anchor, const segment& target) {
        if (num==NPAIRS) { dump(); }

        int awidth=anchor.width;
        int twidth=target.width;           
		if (awidth<0 || twidth<0) { throw std::runtime_error("alignment lengths should be positive"); }
        if (anchor.reverse) { awidth *= -1; } 
        if (target.reverse) { twidth *= -1; }

        ai[num]=anchor.fragid+1; // Get back to 1-indexing.
        ti[num]=target.fragid+1; 
        ap[num]=anchor.pos;
        tp[num]=target.pos;
        al[num]=awidth;
        tl[num]=twidth;
        ++num;
        return;
    }

    void dump() {
        if (!num) { return; }
        if (saved) {
            out=std::fopen(path.c_str(), "a");
        } else {
            out=std::fopen(path.c_str(), "w"); // Overwrite any existing file, just to be safe.
        }
        if (out==NULL) {
            std::stringstream err;
            err << "failed to open output file at '" << path << "'"; 
            throw std::runtime_error(err.str());
        }
        for (size_t i=0; i<num; ++i) {
            fprintf(out, "%i\t%i\t%i\t%i\t%i\t%i\n", ai[i], ti[i], ap[i], tp[i], al[i], tl[i]);
        }
        std::fclose(out);
        num=0;
        saved=true;
        return;
    }

    size_t num;
    const size_t NPAIRS;
    std::deque<int> ai, ti, ap, tp, al, tl;
    std::string path;
    FILE * out; 
    bool saved;
};

/************************
 * Main loop.
 ************************/

SEXP internal_loop (const base_finder * const ffptr, status (*check_self_status)(const segment&, const segment&), const check_invalid_chimera * const icptr,
        SEXP chr_converter, SEXP bamfile, SEXP prefix, SEXP storage, SEXP chimera_strict, SEXP minqual, SEXP do_dedup) {

    // Checking input values.
    Rcpp::String bampath=check_string(bamfile, "BAM file path");
    Rcpp::String oprefix=check_string(prefix, "output prefix");
    const bool rm_invalid=check_logical_scalar(chimera_strict, "chimera removal specification");
    const bool rm_dup=check_logical_scalar(do_dedup, "duplicate removal specification");
    const int minq=check_integer_scalar(minqual ,"minimum mapping quality");
	const bool rm_min=!ISNA(minq);
    const size_t stored_pairs=check_integer_scalar(storage, "number of stored pairs");

    Bamfile input(bampath.get_cstring());
    
    // Initializing the chromosome conversion table (to get from BAM TIDs to chromosome indices in the 'fragments' GRanges).
	const size_t nc=ffptr->nchrs();
    Rcpp::IntegerVector converter(chr_converter);
    const int nbamc=converter.size();
    if (nbamc > int(nc)) { throw std::runtime_error("more chromosomes in the BAM file than in the fragment list"); }
    for (int i=0; i<nbamc; ++i) {
        if (converter[i]==NA_INTEGER || converter[i] < 0 || converter[i] >= int(nc)) { throw std::runtime_error("conversion indices out of range"); }
    }
    
   	// Constructing output containers
	std::vector<std::deque<OutputFile> > collected(nc);
	for (size_t i=0; i<nc; ++i) { 
        for (size_t j=0; j<=i; ++j) { 
            collected[i].push_back(OutputFile(oprefix.get_cstring(), i, j, stored_pairs));
        }
    }
	int single=-1; // First one always reported as a singleton, as qname is empty.
	int total=0, dupped=0, filtered=0, mapped=0;
	int dangling=0, selfie=0;
	int total_chim=0, mapped_chim=0, multi_chim=0, inv_chimeras=0;

    std::string qname="";
    std::deque<segment> read1, read2;
    while (1) {
        bool isempty=true;
        int nsegments=0;
        bool isdup=false;
        bool firstunmap=true, secondunmap=true;
        bool hasfirst=false, hassecond=false;
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
			const bool isfirst=bool((input.read -> core).flag & BAM_FREAD1);
			if (isfirst) { hasfirst=true; }
			else { hassecond=true; }
            
			// Checking how we should proceed; whether we should bother adding it or not.
			const bool curdup=bool((input.read -> core).flag & BAM_FDUP);
			const bool curunmap=(bool((input.read -> core).flag & BAM_FUNMAP) || (rm_min && (input.read -> core).qual < minq));
            int offset, width;
			parse_cigar(input.read, offset, width);
            if (offset==0 && width > 0) { 
                if (curdup) { isdup=true; } // defaults to 'false' unless we have a definitive setting of markingness.
                if (!curunmap) { (isfirst ? firstunmap : secondunmap)=false; } // defaults to 'true' unless we know it's mapped (unmapped reads get width=0 and won't reach here). 
            }

			// Checking which deque to put it in, if we're going to keep it.
            if (! (curdup && rm_dup) && ! curunmap) {
                const int32_t& curtid=(input.read -> core).tid;
                if (curtid==-1 || curtid >= nbamc) {
                    std::stringstream err;
                    err << "tid for read '" << bam_get_qname(input.read) << "' out of range of BAM header";
                    throw std::runtime_error(err.str());
                } 

                segment current(converter[curtid], // Chromosome ID
                                (input.read->core).pos + 1, // Code assumes 1-based index for base position.
                                bool(bam_is_rev(input.read)), // Specifies if reverse.
                                offset, width);

                std::deque<segment>& current_reads=(isfirst ? read1 : read2);
                if (offset==0) { current_reads.push_front(current); } 
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
        const bool ischimera=(nsegments > 2);
		if (ischimera) { ++total_chim; }
		if (isdup) { ++dupped; }
        const bool isunmap=(firstunmap | secondunmap);
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
			current.fragid=ffptr->find_fragment(current);
		}
		for (size_t i2=0; i2<read2.size(); ++i2) {
			segment& current=read2[i2];
			current.fragid=ffptr->find_fragment(current);
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
				if (read1.front().get_5pos() > read2.front().get_5pos()) { // Using the 5' ends to determine ordering.
					anchor=true;
				}
			}
		}
		const segment& anchor_seg=(anchor ? read1.front() : read2.front());
		const segment& target_seg=(anchor ? read2.front() : read1.front());   
        collected[anchor_seg.chrid][target_seg.chrid].add(anchor_seg, target_seg);
	}

    // Dumping any leftovers that are still present.
    for (size_t i=0; i<nc; ++i) { 
        for (size_t j=0; j<=i; ++j) { 
            collected[i][j].dump();
        }
    }

    // Saving all file names.
    Rcpp::List filepaths(nc);
    for (size_t i=0; i<nc; ++i) {
        Rcpp::StringVector outpaths(i+1);
        for (size_t j=0; j<=i; ++j) {
            if (collected[i][j].saved) {
                outpaths[j]=collected[i][j].path;
            }
        } 
        filepaths[i]=outpaths;
    }

    return Rcpp::List::create(filepaths, 
        Rcpp::IntegerVector::create(total, dupped, filtered, mapped), // Dumping mapping diagnostics.
        Rcpp::IntegerVector::create(dangling, selfie), // Dumping self-interaction diagnostics.
        Rcpp::IntegerVector::create(single), // Dumping the number designated 'single'.
        Rcpp::IntegerVector::create(total_chim, mapped_chim, multi_chim, inv_chimeras)); // Dumping chimeric diagnostics.
}

SEXP report_hic_pairs (SEXP start_list, SEXP end_list, SEXP chrconvert, SEXP bamfile, SEXP outfile, SEXP storage, 
        SEXP chimera_strict, SEXP chimera_span, SEXP minqual, SEXP do_dedup) {
    BEGIN_RCPP 
	fragment_finder ff(start_list, end_list);
	
	check_invalid_by_fragid invfrag; // Bit clunky to define both, but easiest to avoid nested try/catch.
	check_invalid_by_dist invdist(chimera_span);
	const check_invalid_chimera* invchim=NULL;
	if (invdist.get_span()==NA_INTEGER) { invchim=&invfrag; } 
	else { invchim=&invdist; }
	
	return internal_loop(&ff, &get_status, invchim, chrconvert, bamfile, outfile, storage, chimera_strict, minqual, do_dedup);
    END_RCPP;
}

/************************
 * Repeated loop for DNase Hi-C.
 ************************/

class simple_finder : public base_finder {
public:
	simple_finder(SEXP);
    size_t nchrs() const;
	int find_fragment(const segment&) const;
private:
    std::vector<int> chrlens;
};

size_t simple_finder::nchrs() const {
    return chrlens.size();
}

simple_finder::simple_finder(SEXP clens) { 
    Rcpp::IntegerVector _clens(clens);
    chrlens.insert(chrlens.end(), _clens.begin(), _clens.end());
	return;	
}

int simple_finder::find_fragment(const segment& current) const {
	if (current.reverse && current.get_5pos() > chrlens[current.chrid]) { 
        Rcpp::warning("read aligned off end of chromosome"); 
    }
    return 0;
}

status no_status_check (const segment& left, const segment& right) {
	/* Fragment IDs have no concept in DNase Hi-C, so automatic
	 * detection of self-circles/dangling ends is impossible.
	 */
	(void)left;
	(void)right; // Just to avoid unused warnings, but maintain compatibility.
	return NEITHER;
}

SEXP report_hic_binned_pairs (SEXP chrlens, SEXP chrconvert, SEXP bamfile, SEXP outfile, SEXP storage,
        SEXP chimera_strict, SEXP chimera_span, SEXP minqual, SEXP do_dedup) {
    BEGIN_RCPP
	simple_finder ff(chrlens);
	check_invalid_by_dist invchim(chimera_span);
	return internal_loop(&ff, &no_status_check, &invchim, chrconvert, bamfile, outfile, storage, chimera_strict, minqual, do_dedup);
    END_RCPP
}

/********************
 * Testing functions.
 *******************/

SEXP test_parse_cigar (SEXP incoming) {
    BEGIN_RCPP

    Rcpp::String bampath=check_string(incoming, "BAM file path");
    Bamfile input(bampath.get_cstring());
    if (sam_read1(input.in, input.header, input.read)<0) { 
        throw std::runtime_error("BAM file is empty");
    } 
   
    Rcpp::IntegerVector output(2);
    parse_cigar(input.read, output[1], output[0]);
    return output;
    END_RCPP
}

SEXP test_fragment_assign(SEXP starts, SEXP ends, SEXP chrs, SEXP pos, SEXP rev, SEXP len) {
    BEGIN_RCPP
	fragment_finder ff(starts, ends);

    Rcpp::IntegerVector _chrs(chrs), _pos(pos), _len(len);
    Rcpp::LogicalVector _rev(rev);
	const int n=_chrs.size();
	if (n!=_pos.size() || n!=_len.size() || n!=_rev.size()) { throw std::runtime_error("length of data vectors are not consistent"); }

    Rcpp::IntegerVector output(n);
	for (int i=0; i<n; ++i) {
        segment current(_chrs[i], _pos[i], bool(_rev[i]), 0, _len[i]); 
		output[i]=ff.find_fragment(current)+1;
	}
	
	return output;
    END_RCPP
}

