cutGenome <- function(bs, pattern, overhang=4L) 
# This finds the target cut sites in the genome. It currently only searches the
# sense strand, which is fine because if the patterns is an inverse palindrome.
# Otherwise, there may be some problems as you'd need to search the reverse
# strand as well.
#
# written by Aaron Lun
# a long time ago. 
# last modified 22 March 2014
{
	if (nchar(pattern)%%2L!=0) { stop("recognition site must be even in size") }
	ps <- DNAString(pattern)
	if (reverseComplement(ps)!=ps) { stop("recognition site must be an inverse palindrome") }
	overhang <- as.integer(overhang)
	if (overhang > nchar(pattern) || overhang < 0L || overhang%%2L!=0) { stop("overhang must be a non-negative even integer that is not greater than pattern length") }
	remainder <- (nchar(pattern)-overhang)/2L

	original <- list()
	if (is(bs, "BSgenome")) {
		ref.names <- seqnames(bs)
		gen <- genome(bs)
	} else {
		bs <- readDNAStringSet(bs)
		ref.names <- names(bs)
		gen <- NA
	}	

	for (chr in ref.names) {
       	x <- matchPattern(pattern, bs[[chr]])
		match.start <- start(x)
		if (is.unsorted(match.start)) { match.start <- sort(match.start) }

		chrlen <- length(bs[[chr]])
   		starts <- c(1L, match.start+remainder)
   		ends <- c(match.start+remainder-1L+overhang, chrlen)
		original[[chr]] <- GRanges(chr, IRanges(starts, ends))
	}
		
	names(original) <- NULL
	suppressWarnings(original <- do.call(c, original))
    seqlevels(original) <- ref.names
    suppressWarnings(seqinfo(original) <- Seqinfo(ref.names, seqlengths=seqlengths(bs), genome=gen))
	return(original)
}

