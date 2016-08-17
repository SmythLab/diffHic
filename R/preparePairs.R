preparePairs <- function(bam, param, file, dedup=TRUE, minq=NA, ichim=TRUE, chim.dist=NA, output.dir=NULL, storage=5000)
# This function prepares Hi-C data by stripping out all valid pairs from the BAM file and
# returning a table describing the interacting fragments of that pair. Diagnnostic data is
# also returned describing various bits and pieces of hiC quality.
#
# written by Aaron Lun
# created 30 May 2013
# last modified 9 December 2015
{
	# Preparing cuts; start positions, end positions, index in 'fragments', segmented by chromosome.
	# Anchor order is defined by the order of chromosomes in 'fragments'; earlier chromosomes
	# are designated as the second anchor when compared to later chromosomes.
	scuts <- ecuts <- list()
	boost.idx<- list()
	fragments <- param$fragments
	frag.data <- .splitByChr(fragments)
	chrs <- frag.data$chr

	curends<-end(fragments)
	curstarts<-start(fragments)
	for (x in seq_along(chrs)) {
		curdex <- frag.data$first[x]:frag.data$last[x]
		scuts[[x]] <- curstarts[curdex]
		ecuts[[x]] <- curends[curdex]
		if (x==1L) {
			boost.idx[[chrs[x]]] <- 0L
		} else {
			boost.idx[[chrs[x]]] <- frag.data$last[x-1L]
		}
	}

	# Checking consistency between SAM chromosome lengths and the ones in the cuts.
	chromosomes<-scanBamHeader(bam)[[1]]$targets
    bam.chrs <- names(chromosomes)
    m <- match(bam.chrs, chrs)
	if (any(is.na(m))) { stop("missing chromosomes in cut site list") }
    for (x in seq_along(bam.chrs)) {
		if (chromosomes[x]!=tail(ecuts[[m[x]]], 1)) {
			stop("length of ", bam.chrs[x], " is not consistent between BAM file and fragment ranges")
		}
	}

	# Enforcing input types.
	minq <- as.integer(minq)
	ichim <- as.logical(ichim)
	chim.dist <- as.integer(chim.dist)
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

    # Calling the C++ code that does everything.
	out <- .Call(cxx_report_hic_pairs, scuts, ecuts, m-1L, path.expand(bam), prefix, storage, !ichim, chim.dist, minq, dedup)
    if (is.character(out)) { stop(out) }
    .process_output(out, file, chrs, boost.idx)
}

.splitByChr <- function(ranges)
# Gets the start and end for each chromosome in the sorted GRanges. 
{
	chrs <- as.character(runValue(seqnames(ranges)))
	if (anyDuplicated(chrs)) { stop("ranges for each chromosome should be consecutive") }
	ref.len <- runLength(seqnames(ranges))
	end.index <- cumsum(ref.len)
	start.index <- end.index - ref.len + 1L
	names(end.index) <- names(start.index) <- chrs
	return(list(chr=chrs, first=start.index, last=end.index))
}

.process_output <- function(c_out, file, chrs, chr.start) {
    .initializeH5(file)
    for (a1.dex in seq_along(c_out[[1]])) { 
        curnames <- c_out[[1]][[a1.dex]]
        not.empty <- curnames!=""
        if (!any(not.empty)) { next }
        anchor1 <- chrs[a1.dex]
        .addGroup(file, anchor1)

		for (a2.dex in which(not.empty)) { 
            anchor2 <- chrs[a2.dex]
            current.file <- curnames[a2.dex]
			out <- read.table(current.file, header=FALSE, colClasses="integer")
			colnames(out) <- c("anchor1.id", "anchor2.id", "anchor1.pos", "anchor2.pos", "anchor1.len", "anchor2.len")

			out$anchor1.id <- out$anchor1.id+chr.start[[anchor1]]
			out$anchor2.id <- out$anchor2.id+chr.start[[anchor2]]
			out <- out[order(out$anchor1.id, out$anchor2.id),,drop=FALSE]
			rownames(out)<-NULL
			.writePairs(out, file, anchor1, anchor2)
		}
	}

    c_out <- c_out[-1]
    names(c_out) <- c("pairs", "same.id", "singles", "chimeras")
	names(c_out$pairs) <-c("total", "marked", "filtered", "mapped")
 	names(c_out$same.id) <- c("dangling", "self.circle")
	names(c_out$chimeras) <- c("total", "mapped", "multi", "invalid")
    return(c_out)
}

