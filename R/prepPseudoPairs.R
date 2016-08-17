prepPseudoPairs <- function(bam, param, file, dedup=TRUE, minq=NA, ichim=TRUE, chim.span=1000, output.dir=NULL, storage=5000L)
# This function acts the same as preparePairs, but it assumes that you're
# putting things into contiguous bins across the genome. The idea is to
# allow DNase-digested Hi-C experiments to fit in the pipeline, where reads
# are assigned to pseudo-restriction fragments.
#
# written by Aaron Lun
# created 27 March 2015
# last modified 9 December 2015
{
	fragments <- param$fragments
	n.per.chr <- runLength(seqnames(fragments))
	frag.data <- .splitByChr(fragments)

	chrs <- frag.data$chr
	last.in.chr <- frag.data$last
	before.first <- as.list(c(0L, last.in.chr[-length(chrs)]))
	names(before.first) <- chrs

	bin.width <- max(width(fragments))
	if (!all(bin.width==width(fragments)[-last.in.chr])) {
		stop("pseudo-fragments should be constant size") 
	}

	# Checking consistency between SAM chromosome lengths and the ones in the cuts.
	chromosomes<-scanBamHeader(bam)[[1]]$targets
    bam.chrs <- names(chromosomes)
    m <- match(bam.chrs, chrs)
    if (any(is.na(m))) { stop("missing chromosomes in cut site list") }
	for (x in seq_along(bam.chrs)) {
		if (chromosomes[x]!=end(fragments)[last.in.chr[m[x]]]) {
			stop("length of ", bam.chrs[x], " is not consistent between BAM file and fragment ranges")
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

    out <- .Call(cxx_report_hic_binned_pairs, n.per.chr, bin.width, m-1L, path.expand(bam), prefix, storage, !ichim, chim.span, minq, dedup)
    if (is.character(out)) { stop(out) }
    final <- .process_output(out, file, chrs, before.first)
    final$same.id <- NULL
    return(final)
}

segmentGenome <- function(bs, size=500) 
# Segments the genome into pseudo-fragments, for use in assigning reads 
# from a DNase Hi-C experiment.
#
# written by Aaron Lun
# created 27 March 2015
{
	if (is(bs, "BSgenome")) {
		ref.len <- seqlengths(bs)
		gen <- genome(bs)
	} else if (is.character(bs)) {
		bs <- readDNAStringSet(bs)
		ref.len <- width(bs)
		names(ref.len) <- names(bs)
		gen <- NA
	} else {
		ref.len <- bs
		gen <- NA
	}

	everything <- list()
	for (chr in names(ref.len)) {
		bin.dex <- seq_len(ceiling(ref.len[[chr]]/size))
		current <- GRanges(chr, IRanges((bin.dex - 1L)*size + 1L, bin.dex*size))
		if (end(current)[length(current)] > ref.len[[chr]]) { 
			end(current)[length(current)] <- ref.len[[chr]]
		}
		everything[[length(everything)+1L]] <- current
	}
	suppressWarnings(everything <- do.call(c, everything))
	seqlengths(everything) <- ref.len
	return(everything)
}
