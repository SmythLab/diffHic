prepPseudoPairs <- function(bam, param, file, dedup=TRUE, minq=NA, ichim=TRUE, chim.span=1000, output.dir=NULL, storage=5000L)
# This function acts the same as preparePairs, but it assumes that you're
# putting things into contiguous bins across the genome. The idea is to
# allow DNase-digested Hi-C experiments to fit in the pipeline, where reads
# are assigned to pseudo-restriction fragments.
#
# written by Aaron Lun
# created 27 March 2015
# last modified 16 March 2017
{
	fragments <- param$fragments
    if (length(fragments)) { 
        stop("fragments should be empty for DNase-C experiments")
    }
    chrlens <- seqlengths(fragments)
    if (length(chrlens)==0L) { 
        stop("seqlengths were not specified in fragments")
    }
    chrs <- names(chrlens)
    before.first <- as.list(integer(length(chrs))-1L) # to undo 1-indexing.
    names(before.first) <- chrs

	# Checking consistency between SAM chromosome lengths and the ones in the cuts.
	chromosomes<-scanBamHeader(bam)[[1]]$targets
    bam.chrs <- names(chromosomes)
    m <- match(bam.chrs, chrs)
    if (any(is.na(m))) { stop("missing chromosomes in cut site list") }
	for (x in seq_along(bam.chrs)) {
		if (chromosomes[x]!=chrlens[m[x]]) { 
			stop("length of ", bam.chrs[x], " is not consistent between BAM file and fragments")
		}
	}

	# Enforcing input types.
	minq <- as.integer(minq)
	ichim <- as.logical(ichim)
	chim.span <- as.integer(chim.span)
	dedup <- as.logical(dedup)
    storage <- as.integer(storage)
    if (storage <= 1L) { 
        stop("'storage' must be a positive integer")
    }

    # Setting up the output directory.
    if (is.null(output.dir)) { 
        output.dir <- tempfile(tmpdir=".")
    } else {
        output.dir <- path.expand(output.dir)
    }
    if (!dir.create(output.dir)) {
        stop("failed to create output directory")
    }
    on.exit({ unlink(output.dir, recursive=TRUE) })
    prefix <- file.path(output.dir, "")

    out <- .Call(cxx_report_hic_binned_pairs, chrlens, m-1L, path.expand(bam), prefix, storage, !ichim, chim.span, minq, dedup)
    if (is.character(out)) { stop(out) }
    final <- .process_output(out, file, chrs, before.first)
    final$same.id <- NULL
    return(final)
}

segmentGenome <- function(bs) 
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
    
    nothing <- GRanges()
    seqlevels(nothing) <- names(ref.len)
	seqlengths(nothing) <- ref.len
	return(nothing)
}
