domainDirections <- function(files, param, width=50000, span=10)
# This function computes the directionality index for each bin in the genome,
# in order to find domain boundaries.
#
# written by Aaron Lun
# created 27 May 2015
{
	nlibs <- length(files)
	if (nlibs==0) { 
		stop("number of libraries must be positive")
	} else if (width < 0) { 
		stop("width must be a non-negative integer")
	} 
	width<-as.integer(width) 
	span <- as.integer(span) 
	fragments <- param$fragments
	new.pts <- .getBinID(fragments, width)
								
	# Setting up ranges for the fragments and bins.
	chrs <- seqlevelsInUse(fragments)
	frag.by.chr <- .splitByChr(fragments)
	bin.by.chr <- .splitByChr(new.pts$region)
	
	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap
	
	# Running through each pair of chromosomes.
    nlibs <- length(files)
	upcount <- downcount <- matrix(0L, length(new.pts$region), nlibs)

	overall <- .loadIndices(files, chrs, restrict)
	for (chr in names(overall)) {
		current <- overall[[chr]]
		if (!(chr %in% names(current))) { next }

		pairs <- .baseHiCParser(current[[chr]], files, chr, chr,
			chr.limits=frag.by.chr, discard=discard, cap=cap)
		first.index <- bin.by.chr$first[[chr]]
		last.index <- bin.by.chr$last[[chr]]
	
		out <- .Call(cxx_directionality, pairs, new.pts$id, span, first.index, last.index)
		if (is.character(out)) { stop(out) }
		if (!length(out[[1]])) { next }

        pnts <- first.index:last.index
		downcount[pnts,] <- out[[1]]
		upcount[pnts,] <- out[[2]]
	}

    # Return an RSE with up and down counts.
    # No total counts, because we don't load every chromosome pair - might as well call totalCounts() externally if required.
    return(SummarizedExperiment(SimpleList(up=upcount, down=downcount), 
        new.pts$region, metadata=list(param=param, span=span, width=width)))                                      
}

