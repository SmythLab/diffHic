squareCounts <- function(files, param, width=50000, filter=1L)
# This function collates counts across multiple experiments to get the full set of results. It takes 
# a list of lists of lists of integer matrices (i.e. a list of the outputs of convertToInteractions) and
# then compiles the counts into a list object for output. 
#
# written by Aaron Lun
# some time ago
# last modified 29 April 2015
{
	nlibs <- length(files)
	if (nlibs==0) { 
		stop("number of libraries must be positive")
	} else if (width < 0) { 
		stop("width must be a non-negative integer")
	} 
	width <- as.integer(width) 
	filter <- as.integer(filter) 
	fragments <- param$fragments
	new.pts <- .getBinID(fragments, width)
	
	# Setting up ranges for the fragments and bins.
	chrs <- seqlevels(fragments)
	frag.by.chr <- .splitByChr(fragments)
	bin.by.chr <- .splitByChr(new.pts$region)
		
	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Output vectors.
	full.sizes <- integer(nlibs)
	out.counts <- list(matrix(0L, 0, nlibs))
	out.a <- out.t <- list(integer(0))
	idex <- 1L

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
    for (anchor in names(overall)) {
        current <- overall[[anchor]]
		for (target in names(current)) {

			# Extracting counts and checking them.
			pairs <- .baseHiCParser(current[[target]], files, anchor, target, 
				chr.limits=frag.by.chr, discard=discard, cap=cap)
			full.sizes <- full.sizes + sapply(pairs, FUN=nrow)
			
			# Aggregating them in C++ to obtain count combinations for each bin pair.
			out <- .Call(cxx_count_patch, pairs, new.pts$id, filter, 
				bin.by.chr$first[[target]], bin.by.chr$last[[target]])
			if (is.character(out)) { stop(out) }
			if (!length(out[[1]])) { next }

			# Storing counts and locations. 
			if (any(out[[1]] < out[[2]])) { stop("anchor ID should not be less than target ID") }
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

	return(DIList(counts=out.counts, totals=full.sizes, 
		anchors=out.a, targets=out.t, regions=new.pts$region,
		exptData=List(param=param, width=width)))
}

## PROOF:
# Recall the enforcement of anchor >= target. Bin pairs should technically be
# reflected around the diagonal, to ensure that all points are counted, e.g.,
# if a bin pair overlaps a point in (target, anchor) form but not in (anchor,
# target) form. However, this is not required, as shown below:
# 
# Consider a point (x, y) past the diagonal (i.e., no enforcement), where y >
# x. Assume that this point is covered by our bin pair. This implies that the
# target range of our bin pair `[te, ts]` includes 'y', and the anchor range
# `[ae, as]` includes 'x'. The anchor range must be above the target range
# (i.e., as >= ts, ae >= te). If ts <= y <= te and as <= x <= ae, then you can
# fiddle with this to obtain ts <= x <= te and as <= y <= ae (as x < y), i.e.,
# the point (y, x) is also covered. So, we're guaranteed to have already
# counted anything past the diagonal, meaning that reflection is not needed.

####################################################################################################

.getBinID <- function(fragments, width) 
# Determines which bin each restriction fragment is in. Also records the rounded
# start and stop site for each bin. Returns a set of bin ids for each restriction
# fragment on each chromosome, as well as the coordinates of each bin.
{
	width<-as.integer(width)
	out.ids<-integer(length(fragments))
	out.ranges<-list()
	last<-0L
	frag.data <- .splitByChr(fragments)
	nfrags <- list() 
	
	for (x in 1:length(frag.data$chr)) {
		curindex <- frag.data$first[x]:frag.data$last[x]
		curf <- fragments[curindex]
		mids <- (start(curf)+end(curf))/2
		bin.id <- as.integer((mids-0.1)/width)+1L 
		# The '-0.1' in the preceding step reduces 'mids' that are exact multiples 
		# of 'width', so each bin is from (n*width, (n+1)*width] for integer 'n'.

		processed <- rle(bin.id)
		ns <- length(processed$value)
		processed$values <- 1:ns
		nfrags[[x]] <- processed$length
		out.ids[curindex] <- inverse.rle(processed)+last
		
		endx <- cumsum(processed$length)
		startx <- rep(1L, ns)
		if (ns>=2L) { startx[-1] <- endx[-ns]+1L }
		out.ranges[[x]] <- GRanges(frag.data$chr[x], IRanges(start(curf[startx]), end(curf[endx])))
		last <- last+ns
	}

	# Wrapping up.
	suppressWarnings(out.ranges <- do.call(c, out.ranges))
	seqlevels(out.ranges) <- seqlevels(fragments)
	seqlengths(out.ranges) <- seqlengths(fragments)
	out.ranges$nfrags <- unlist(nfrags)
	return(list(id=out.ids, region=out.ranges))
}

####################################################################################################

.baseHiCParser <- function(ok, files, anchor, target, chr.limits, discard, cap)
# A courtesy function, to assist with loading counts in this function and others.
{
	overall<-list()
	adisc <- discard[[anchor]]
	tdisc <- discard[[target]]
	do.cap <- !is.na(cap)

	for (x in 1:length(ok)) {
		if (!ok[x]) { 
			overall[[x]] <- data.frame(anchor.id=integer(0), target.id=integer(0))
		} else {
			out <- .getPairs(files[x], anchor, target)
	
			# Checking fidelity, figuring out which ones to throw out.
			check <- .Call(cxx_check_input, out$anchor.id, out$target.id)
			if (is.character(check)) { stop(check) }

			# Checking that we're all on the right chromosome.
			if (nrow(out)) { 
				if (max(out$anchor.id) > chr.limits$last[[anchor]] || 
						min(out$anchor.id) < chr.limits$first[[anchor]]) { 
					stop("anchor index outside range of fragment object") 
				}
				if (max(out$target.id) > chr.limits$last[[target]] || 
						min(out$target.id) < chr.limits$first[[target]]) { 
					stop("target index outside range of fragment object") 
				}
			}

			# Overlapping with those in the discard intervals.
			if (!is.null(adisc) || !is.null(tdisc)) {
				a.hits <- t.hits <- FALSE
 			    if (!is.null(adisc)) {
					a.hits <- overlapsAny(IRanges(out$anchor.pos, out$anchor.pos+abs(out$anchor.len)-1L), adisc, type="within")
				}
				if (!is.null(tdisc)) { 
					t.hits <- overlapsAny(IRanges(out$target.pos, out$target.pos+abs(out$target.len)-1L), tdisc, type="within")
				}
				out <- out[!a.hits & !t.hits,,drop=FALSE]
			}

			# Removing read pairs above the cap for each restriction fragment pair.
			if (do.cap) { 
				capped <- .Call(cxx_cap_input, out$anchor.id, out$target.id, cap)
				if (is.character(capped)) { stop(capped) }
				out <- out[capped,]
			}

			dim(out$anchor.id) <- dim(out$target.id) <- NULL
			overall[[x]] <- out[,c("anchor.id", "target.id")]
		}
	}
	return(overall)
}

####################################################################################################

.splitDiscards <- function(discard) 
# Splits the discard GRanges into a list of constituent chromosomes,
# along with IRanges for everything. This allows easy access tot he
{
	if (is.null(discard) || length(discard)==0) { return(NULL) }
	discard <- sort(discard)
	all.chrs <- as.character(runValue(seqnames(discard)))
	all.len <- runLength(seqnames(discard))
	chr.ends <- cumsum(all.len)
	chr.starts <- c(1L, chr.ends[-length(chr.ends)]+1L)

	output <- list()
	for (i in 1:length(all.chrs)) {
		chr <- all.chrs[i]
		ix <- chr.starts[i]:chr.ends[i]
		output[[chr]] <- reduce(ranges(discard[ix]))
	}

	return(output)
}

