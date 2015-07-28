neighborCounts <- function(files, param, width=50000, filter=1L, flank=NULL, exclude=NULL, prior.count=NULL)
# This does the same thing as squareCounts, except that it simultaneously computes the 
# filter statistic for each extracted bin pair. This has lower memory requirements as
# it doesn't need to hold the entire `filter=1` output in memory at once.
#
# written by Aaron Lun
# created 21 May 2015
# last modified 22 July 2015
{
	nlibs <- length(files)
	if (nlibs==0L) {
		stop("number of libraries must be positive")
	} else if (width < 0) { 
		stop("width must be a non-negative integer")
	} 
	width<-as.integer(width) 
	filter<-as.integer(filter) 
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
	out.counts <- list(matrix(0L, 0, nlibs))
	out.a <- out.t <- list(integer(0))
	out.filter <- list(double(0))
	idex <- 1L
	
	# Getting the full library sizes and computing offsets.
	full.sizes <- totalCounts(files, param)
	offsets <- log(full.sizes) - log(mean(full.sizes))
	maxit <- as.integer(formals(mglmOneGroup)$maxit)
	tol <- formals(mglmOneGroup)$tol
	disp <- 0.05

	# Other stuff related to calculation of the neighborhood regions.	
	if (is.null(flank)) { flank <- formals(enrichedPairs)$flank }
	if (is.null(exclude)) { exclude <- formals(enrichedPairs)$exclude }
	flank <- as.integer(flank)
	exclude <- as.integer(exclude)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	if (exclude < 0L) { stop("exclude width must be a positive integer") }
	if (flank <= exclude) { stop("exclude width must be less than the flank width") }
	if (is.null(prior.count)) { prior.count <- formals(enrichedPairs)$prior.count } 
	prior.count <- as.double(prior.count)    

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
    for (anchor in names(overall)) {
        current <- overall[[anchor]]
		for (target in names(current)) {
			pairs <- .baseHiCParser(current[[target]], files, anchor, target, 
				chr.limits=frag.by.chr, discard=discard, cap=cap)
			
			# Aggregating counts in C++ to obtain count combinations for each bin pair.
			out <- .Call(cxx_count_background, pairs, new.pts$id, flank, exclude, filter, 
				bin.by.chr$first[[target]], bin.by.chr$last[[target]], bin.by.chr$first[[anchor]], bin.by.chr$last[[anchor]],
				maxit, tol, offsets, disp, prior.count)
			if (is.character(out)) { stop(out) }
			if (!length(out[[1]])) { next }

			# Storing counts and locations. 
			if (any(out[[1]] < out[[2]])) { stop("anchor ID should not be less than target ID") }
			out.a[[idex]] <- out[[1]]
 			out.t[[idex]] <- out[[2]]
			out.counts[[idex]] <- out[[3]]
			out.filter[[idex]] <- out[[4]]
			idex<-idex+1L
		}
	}

	# Collating all the other results.
	out.a <- unlist(out.a)
	out.t <- unlist(out.t)
	out.f <- log2(unlist(out.filter))
	out.counts <- do.call(rbind, out.counts)

	return(list(interaction=DIList(counts=out.counts, totals=full.sizes, 
		anchors=out.a, targets=out.t, regions=new.pts$region,
		exptData=List(param=param, width=width)), enrichment=out.f))
}


