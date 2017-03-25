prepPseudoPairs <- function(bam, param, file, dedup=TRUE, minq=NA, ichim=TRUE, chim.span=1000, output.dir=NULL, storage=5000L)
# This function acts the same as preparePairs, but it assumes that you're
# putting things into contiguous bins across the genome. The idea is to
# allow DNase-digested Hi-C experiments to fit in the pipeline, where reads
# are assigned to pseudo-restriction fragments.
#
# written by Aaron Lun
# created 27 March 2015
# last modified 20 March 2017
{
    .Deprecated("preparePairs")
    preparePairs(bam, param, file, dedup=dedup, minq=minq, ichim=ichim, chim.dist=chim.span, output.dir=output.dir, storage=storage)
}

segmentGenome <- function(bs) {
    .Deprecated("emptyGenome")
    emptyGenome(bs)
}

emptyGenome <- function(bs) 
# Returns an empty 'fragments' but with filled-out seqlengths,
# for use in assigning reads from a DNase Hi-C experiment.
#
# written by Aaron Lun
# created 27 March 2015
# last modified 17 March 2017
{
	if (is(bs, "BSgenome")) {
		ref.len <- seqlengths(bs)
	} else if (is.character(bs)) {
		bs <- readDNAStringSet(bs)
		ref.len <- width(bs)
		names(ref.len) <- names(bs)
	} else {
		ref.len <- bs
	}
    
    GRanges(seqlengths=ref.len)
}
