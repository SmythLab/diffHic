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
	chrs <- seqlevels(fragments)
	frag.by.chr <- .splitByChr(fragments)
	bin.by.chr <- .splitByChr(new.pts$region)
	
	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap
	
	# Getting the full library sizes and computing offsets.
	full.sizes <- totalCounts(files, param)
	offsets <- log(full.sizes) - log(mean(full.sizes))
	maxit <- as.integer(formals(mglmOneGroup)$maxit)
	tol <- formals(mglmOneGroup)$tol
	disp <- 0.05
	
	# Running through each pair of chromosomes.
	upcount <- downcount <- numeric(length(new.pts$region))
	overall <- .loadIndices(files, chrs, restrict)
	for (chr in names(overall)) {
		current <- overall[[chr]]
		if (!(chr %in% names(current))) { next }

		pairs <- .baseHiCParser(current[[chr]], files, chr, chr,
			chr.limits=frag.by.chr, discard=discard, cap=cap)
		first.index <- bin.by.chr$first[[chr]]
		last.index <- bin.by.chr$last[[chr]]
	
		out <- .Call(cxx_directionality, pairs, new.pts$id, span, 
			first.index, last.index, maxit, tol, offsets, disp)
		if (is.character(out)) { stop(out) }
		if (!length(out[[1]])) { next }

		upcount[first.index:last.index] <- out[[2]]
		downcount[first.index:last.index] <- out[[1]]
	}

	# Computing the significance of deviations.
	new.pts$region$Up <- upcount
	new.pts$region$Down <- downcount
	return(new.pts$region)
}

