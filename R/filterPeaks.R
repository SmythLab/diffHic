enrichedPairs <- function(data, flank=5, exclude=0, prior.count=2, abundances=NULL)
# This function identifies the highest-abundance neighbor in the interaction space
# for each bin pair in `data`. The aim is to compare the abundance of each element
# with the abundance of its neighbor. 
#
# written by Aaron Lun
# created 23 April 2014
# last modified 22 July 2015
{
	flank <- as.integer(flank)
	exclude <- as.integer(exclude)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	if (exclude < 0L) { stop("exclude width must be a positive integer") }
	if (flank <= exclude) { stop("exclude width must be less than the flank width") }
	
	rdata <- .splitByChr(regions(data))
	last.id <- rdata$last
	first.id <- rdata$first
	if (is.null(abundances)) { abundances <- aveLogCPM(asDGEList(data), prior.count=0) }

	# Rescaling to count-level data with at least 6 dp, for stable calculations with integers.
	# This should be enough precision while avoiding overrun of the integer type.
	scaling <- log2(mean(data$totals)/1e6 * ncol(data))
	MULT <- 1e6
	back2count <- function(ab) { 
		newval <- 2^(ab + scaling)
		lower <- as.integer(newval)
		list(val=newval, int=lower, dec=as.integer(round((newval-lower)*MULT)))
	}
	prior.count <- prior.count * ncol(data)	
	
	# Running through each pair of chromosomes.
	np <- nrow(data)
	all.chrs <- as.character(seqnames(regions(data)))
	aid <- anchors(data, id=TRUE)
	by.chr <- split(seq_len(np), all.chrs[aid])
	tid <- targets(data, id=TRUE)
	output <- numeric(np)

	for (anchor in names(by.chr)) {
		next.chr <- by.chr[[anchor]]
		next.chr <- split(next.chr, all.chrs[tid[next.chr]])
		a.len <- last.id[[anchor]] - first.id[[anchor]] + 1L

		for (target in names(next.chr)) {
			current.pair <- next.chr[[target]]
			all.a <- aid[current.pair] - first.id[[anchor]] 
			all.t <- tid[current.pair] - first.id[[target]]
			t.len <- last.id[[target]] - first.id[[target]] + 1L

			rel.ab <- abundances[current.pair]
			converted <- back2count(rel.ab) 

			# Using the quadrant with the maximum average.
			o <- order(all.a, all.t)
			collected <- .Call(cxx_quadrant_bg, all.a[o], all.t[o], 
				converted$int[o], converted$dec[o], MULT, 
				flank, exclude, a.len, t.len, anchor==target)
			if (is.character(collected)) { stop(collected) }
			collected[o] <- collected

			output[current.pair] <- log2((converted$val+prior.count)/(collected+prior.count))
		}
	}

	# Returning the collected log-fold-changes.
	return(output)
}

filterPeaks <- function(data, enrichment, min.enrich=log2(1.5), min.count=5, min.diag=2L, ...)
# This is a wrapper function that takes the enrichment values from filterPeaks and 
# identifies those bin pairs that satisfy the enrichment threshold, and other cut-offs.
# 
# written by Aaron Lun
# created 23 March 2015
# last modified 24 March 2015
{
	keep <- enrichment > min.enrich 
	if (!is.null(min.count)) { 
		ab <- aveLogCPM(asDGEList(data), ...)
		keep <- keep & ab > aveLogCPM(min.count, lib.size=mean(data$totals), ...)
	} 
	if (!is.null(min.diag)) { 
		keep <- keep & (is.na(getDistance(data))  |
			anchors(data, id=TRUE) - targets(data, id=TRUE) >= min.diag)
	}
	return(keep)
}
