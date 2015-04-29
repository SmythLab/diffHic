preparePairs <- function(bam, param, file, dedup=TRUE, yield=1e7, ichim=TRUE, minq=NA)
# This function prepares Hi-C data by stripping out all valid pairs from the BAM file and
# returning a table describing the interacting fragments of that pair. Diagnnostic data is
# also returned describing various bits and pieces of hiC quality.
#
# written by Aaron Lun
# created 30 May 2013
# last modified 29 April 2015
{
	# Preparing cuts; start positions, end positions, index in 'fragments', segmented by chromosome.
	# Anchor/target order is defined by the order of chromosomes in 'fragments'; earlier chromosomes
	# are designated as targets when compared to later chromosomes.
	scuts <- ecuts <- list()
	boost.idx<- list()
	fragments <- param$fragments
	frag.data <- .splitByChr(fragments)
	chrs <- frag.data$chr

	curends<-end(fragments)
	curstarts<-start(fragments)
	for (x in 1:length(chrs)) {
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
	if (!all(names(chromosomes) %in% chrs)) { stop("missing chromosomes in cut site list") }
	for (x in 1:length(chrs)) {
		if (chromosomes[[chrs[x]]]!=tail(ecuts[[x]], 1)) {
			stop("length of ", chrs[x], " is not consistent between BAM file and fragment ranges")
		}
	}

	# Enforcing input types.
	minq <- as.integer(minq)
	ichim <- as.logical(ichim)
	dedup <- as.logical(dedup)

	# Setting up the calling function.
	FUN <- function(read.pair.len, cur.chrs, out) { 
		.Call(cxx_report_hic_pairs, scuts, ecuts, read.pair.len, cur.chrs, out$pos,
			out$flag, out$cigar, out$mapq, !ichim, minq, dedup)
	}

	# Returning collated results.
	.innerPrepLoop(bam=bam, file=file, chrs=chrs, chr.start=boost.idx, FUN=FUN, yield=yield)
}

.splitByChr <- function(ranges)
# Gets the start and end for each chromosome in the sorted GRanges. 
{
	ref.chrs <- as.character(runValue(seqnames(ranges)))
	if (anyDuplicated(ref.chrs)) { stop("ranges for each chromosome should be consecutive") }
	ref.len <- runLength(seqnames(ranges))
	end.index <- cumsum(ref.len)
	start.index <- end.index - ref.len + 1L
	names(end.index) <- names(start.index) <- ref.chrs
	return(list(chr=ref.chrs, first=start.index, last=end.index))
}

.innerPrepLoop <- function(bam, file, chrs, chr.start, FUN, yield)
# This is the function that does the heavy lifting of looping through the file.
# The collation function is specified in FUN. The idea is to allow the prepBinPairs
# function to use the same loop without having to repeat the code.
{
	# Setting up storage vectors for diagnostics and other output.
	diagnostics<- 0L
	same.id <- 0L
	singletons <- 0L
	chimeras <- 0L
	allfiles<-list()
	file.count<-1L

	# Running through all pairs. Note that, despite the yieldSize, elements are extracted so
	# that runs of 'QNAME's are not broken. See:
	# 	https://stat.ethz.ch/pipermail/bioconductor/2013-March/051490.html for more details.
	bf<-open(BamFile(bam, yieldSize=yield, obeyQname=TRUE, index=character(0)))
	dir <- tempfile(tmpdir=".")
	on.exit(unlink(dir, recursive=TRUE))
	dir.create(dir)

	while (1) {
		out <- scanBam(bf, param=ScanBamParam(what=c("qname", "flag", "rname", "pos", "mapq", "cigar")))[[1]]
		if (!length(out[[1]])) { break; }

		# Converting chromosome ids, using the full set of chromosomes (not just those in this yield cycle).
		stopifnot(is.factor(out$rname))
		rematched <- match(levels(out$rname), chrs)
		if (any(is.na(rematched))) { stop("unrecognised chromosomes in the BAM file") }
		cur.chrs <- rematched[as.integer(out$rname)] - 1L

		# Collating read names. We do it here so that any optimization of string comparisons applies
		# here, rather than having to modify the C++ code from strcmp to pointer comparisons.
		read.pair.len <- rle(out$qname)$length
		collated <- FUN(read.pair.len, cur.chrs, out)
		if (is.character(collated)) { stop(collated) }

		# Adding the statistics.
		diagnostics <- diagnostics + collated[[3]]
		same.id <- same.id + collated[[4]]
 		singletons <- singletons + collated[[5]]
		chimeras <- chimeras + collated[[6]]
		out <- NULL

		# Dumping to a temporary place.
		nonempty <- collated[[1]]
		if (!nrow(nonempty)) { next }
		pairdata <- collated[[2]]

		for (i in 1:nrow(nonempty)) {
			anchor <- chrs[nonempty[i,1]]
			target <- chrs[nonempty[i,2]]
			if (is.null(allfiles[[anchor]])) { allfiles[[anchor]] <- list() }
	
			current.file <- allfiles[[anchor]][[target]]
			if (is.null(current.file)) {
				current.file <- file.path(dir, paste0(file.count, ".gz"))
				allfiles[[anchor]][[target]] <- current.file
				file.count <- file.count+1L
			}
	
			fout <- gzfile(current.file, open="ab")
			write.table(file=fout, pairdata[[i]], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
			close(fout)
		}
		collated <- NULL
	}
	close(bf)

	# Parsing the directory, pulling out data.frames. We adjust the a/t indices to get the
	# full values, and we sort them by anchor/target. We then save them into a HDF5 file.
	.initializeH5(file)
	for (anchor in names(allfiles)) {
		tfiles <- allfiles[[anchor]]
		.addGroup(file, anchor)
		for (target in names(tfiles)) {
			current.file <- tfiles[[target]]
			out <- read.table(current.file, header=FALSE, colClasses="integer")
			colnames(out) <- c("anchor.id", "target.id", "anchor.pos", "target.pos", "anchor.len", "target.len")
			out$anchor.id <- out$anchor.id+chr.start[[anchor]]
			out$target.id <- out$target.id+chr.start[[target]]
			out <- out[order(out$anchor.id, out$target.id),,drop=FALSE]
			rownames(out)<-NULL
			.writePairs(out, file, anchor, target)
		}
	}

	# Returning a list of diagnostic components.
	names(diagnostics)<-c("total", "marked", "filtered", "mapped")
 	names(same.id) <- c("dangling", "self.circle")
	names(chimeras) <- c("total", "mapped", "multi", "invalid")
	return(list(pairs=diagnostics,
				same.id=same.id,
				singles=singletons,
				chimeras=chimeras))
}

