marginCounts <- function(files, param, width=50000)
# Gets the marginal counts i.e. sum of counts for each bin or region.
# This is useful to determine the `genomic coverage' of each region,
# based on the number of Hi-C read pairs involving that region.
#
# written by Aaron Lun
# Some time ago.
# last modified 20 March 2015
{
	nlibs <- length(files)
	width <- as.integer(width)
	fragments <- param$fragments

	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

    if (width < 0) { stop("width must be a non-negative integer") }
    new.pts <- .getBinID(fragments, width)
	total.bins <- length(new.pts$region)
	stopifnot(max(new.pts$id)==total.bins) 

    total.counts <- matrix(0L, length(new.pts$region), nlibs)
	full.sizes <- integer(nlibs)
	chrs <- seqlevels(fragments)

    # Running through each pair of chromosomes.
    overall <- .loadIndices(files, chrs, restrict)
    for (anchor in names(overall)) {
		current <- overall[[anchor]]
		for (target in names(current)) {
			if (!.checkIfPairOK(restrict, anchor, target)) { next } 
    
      		pairs <- .baseHiCParser(current[[target]], files, anchor, target, discard=discard, cap=cap)
           	full.sizes <- full.sizes + sapply(pairs, FUN=nrow)

			# Aggregating them in C++ to get the count combinations and location of each bin.
            out <- .Call(cxx_count_marginals, pairs, new.pts$id, total.bins)
            if (is.character(out)) { stop(out) }
			total.counts <- total.counts+out
		}
	}
	
	# Aggregating all elements.
	retained <- which(rowSums(total.counts)>0.5)
	return(DIList(counts=total.counts[retained,,drop=FALSE], totals=full.sizes, 
			anchors=retained, targets=retained, regions=new.pts$region, 
			exptData=List(param=param)))
}

