squareCounts <- function(files, param, width=50000, filter=1L)
# This function collates counts across multiple experiments to get the full set of results. It takes 
# a list of lists of lists of integer matrices (i.e. a list of the outputs of convertToInteractions) and
# then compiles the counts into a list object for output. 
#
# written by Aaron Lun
# some time ago
# last modified 17 March 2017
{
	nlibs <- length(files)
	if (nlibs==0L) {
		stop("number of libraries must be positive")
	} else if (width < 0) { 
		stop("width must be a non-negative integer")
	} 
	width <- as.integer(width) 
	filter <- as.integer(filter) 

    # Setting up the bins.

    # Setting up the other statistics.
    parsed <- .parseParam(param, width=width, bin=TRUE)
    chrs <- parsed$chrs
    frag.by.chr <- parsed$frag.by.chr
    cap <- parsed$cap
    bwidth <- parsed$bwidth
    discard <- parsed$discard
	bin.id <- parsed$bin.id
    bin.region <- parsed$bin.region
	bin.by.chr <- parsed$bin.by.chr
    restrict <- parsed$restrict

	# Output vectors.
	full.sizes <- integer(nlibs)
	out.counts <- list(matrix(0L, 0, nlibs))
	out.a <- out.t <- list(integer(0))
	idex <- 1L

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
    for (anchor1 in names(overall)) {
        current <- overall[[anchor1]]
		for (anchor2 in names(current)) {

			# Extracting counts and checking them.
			pairs <- .baseHiCParser(current[[anchor2]], files, anchor1, anchor2, 
				chr.limits=frag.by.chr, discard=discard, cap=cap, width=bwidth)
			full.sizes <- full.sizes + sapply(pairs, FUN=nrow)
			
			# Aggregating them in C++ to obtain count combinations for each bin pair.
			out <- .Call(cxx_count_patch, pairs, bin.id, filter, 
				bin.by.chr$first[[anchor2]], bin.by.chr$last[[anchor2]])
			if (is.character(out)) { stop(out) }
			if (!length(out[[1]])) { next }

			# Storing counts and locations. 
			if (any(out[[1]] < out[[2]])) { stop("anchor1 ID should not be less than anchor2 ID") }
			out.a[[idex]] <- out[[1]]
 			out.t[[idex]] <- out[[2]]
			out.counts[[idex]] <- out[[3]]
			idex<-idex+1L
		}
	}

	# Collating all the other results.
	out.a <- unlist(out.a)
	out.t <- unlist(out.t)
	out.counts <- do.call(rbind, out.counts)

	return(InteractionSet(list(counts=out.counts), colData=DataFrame(totals=full.sizes), 
		interactions=GInteractions(anchor1=out.a, anchor2=out.t, regions=bin.region, mode="reverse"), 
        metadata=List(param=param, width=width)))
}

## PROOF:
# Recall the enforcement of anchor1 >= anchor2. Bin pairs could technically be
# reflected around the diagonal, to ensure that all points are counted, e.g.,
# if a bin pair overlaps a point in (anchor2, anchor1) form but not in (anchor1,
# anchor2) form. However, this is not required, as shown below:
# 
# Consider a point (x, y) past the diagonal (i.e., no enforcement), where y >
# x. Assume that this point is covered by our bin pair. This implies that the
# anchor2 range of our bin pair `[te, ts]` includes 'y', and the anchor1 range
# `[ae, as]` includes 'x'. The anchor1 range must be above the anchor2 range
# (i.e., as >= ts, ae >= te). If ts <= y <= te and as <= x <= ae, then you can
# fiddle with this to obtain ts <= x <= te and as <= y <= ae (as x < y), i.e.,
# the point (y, x) is also covered. So, we're guaranteed to have already
# counted anything past the diagonal, meaning that reflection is not needed.

