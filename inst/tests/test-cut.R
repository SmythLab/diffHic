###################################################################################################
# This script is designed to test the restriction site finder. We use the E.coli genome to do so 
# because it's short and we won't be spending ages putting it together.	

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

findRestrictionSites <- function(bs, pattern, ref=NULL) {
	require(Biostrings)
	if (nchar(pattern)%%2L!=0L) { stop("Recognition site must be even in size."); }
	ps<-DNAString(pattern);
	if (reverseComplement(ps)!=ps) { stop("Recognition site must be an inverse palindrome."); }
	remainder<-as.integer(nchar(pattern)/2+0.5);
	
	require(BSgenome)
	out<-list();
	for (chr in seqnames(bs)) {
		x<-matchPattern(pattern, bs[[chr]]);
		cuts<-start(x)+remainder-1L;
		out[[chr]]<-c(cuts, length(bs[[chr]]));
	}
	return(out);
}

###################################################################################################

require(BSgenome.Ecoli.NCBI.20080805)

comp<-function(target, overhang=4L) {
	current<-cutGenome(Ecoli, target, overhang=overhang)
	theoretical<-findRestrictionSites(Ecoli, target)

	# Odds and ends.
	overhang<-as.integer(overhang+0.5)
	half.ohang<-as.integer(overhang/2L+0.5)
	half.pattern<-as.integer(nchar(target)/2L)

	# Checking reprocessed cut sites.
	stopifnot(identical(sort(names(theoretical)), sort(seqlevels(current))))
	for (x in names(theoretical)) {
		n<-length(theoretical[[x]])
 	    cuts1<-end(current[x==seqnames(current)])
		cuts1[-n]<-cuts1[-n]-half.ohang
		stopifnot(identical(cuts1, theoretical[[x]]))
	}
	return(head(current))
}

####################################################################################################
# Restriction site must have a 5' overhang.

comp("GGATCC", 4L) # BamHI

comp("GAATTC", 4L) # EcoRI

comp("AACGTT", 2L) # AclI

comp("GCGC", 2L) # HinP1I

comp("GATC", 4L) # MboI

comp("GCGGCCGC", 4L) # NotI

comp("GGCGCGCC", 4L) # AscI

###################################################################################################
# End.
