neighborCounts <- function(files, param, width=50000, filter=1L, flank=NULL, exclude=NULL)
# This does the same thing as squareCounts, except that it simultaneously computes the 
# filter statistic for each extracted bin pair. This has lower memory requirements as
# it doesn't need to hold the entire `filter=1` output in memory at once.
#
# written by Aaron Lun
# created 21 May 2015
# last modified 2 March 2017
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
	chrs <- seqlevelsInUse(fragments)
	frag.by.chr <- .splitByChr(fragments)
	bin.by.chr <- .splitByChr(new.pts$region)
		
	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Output vectors.
	out.counts <- list(matrix(0L, 0, nlibs))
	out.a <- out.t <- list(integer(0))
    full.sizes <- integer(nlibs)

    modes <- .neighbor_locales()
    neighbor.counts <- lapply(modes, FUN=function(x) list(matrix(0L, 0, nlibs)))
    neighbor.N <- lapply(modes, FUN=function(x) list(double(0)))
	idex <- 1L
	
	# Other stuff related to calculation of the neighborhood regions.	
	if (is.null(flank)) { flank <- formals(enrichedPairs)$flank }
	if (is.null(exclude)) { exclude <- formals(enrichedPairs)$exclude }
	flank <- as.integer(flank)
	exclude <- as.integer(exclude)
	if (flank <= 0L) { stop("flank width must be a positive integer") }
	if (exclude < 0L) { stop("exclude width must be a positive integer") }
	if (flank <= exclude) { stop("exclude width must be less than the flank width") }

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
    for (anchor1 in names(overall)) {
        current <- overall[[anchor1]]
		for (anchor2 in names(current)) {
			pairs <- .baseHiCParser(current[[anchor2]], files, anchor1, anchor2, 
				chr.limits=frag.by.chr, discard=discard, cap=cap)
            full.sizes <- full.sizes + sapply(pairs, FUN=nrow)

			# Aggregating counts in C++ to obtain count combinations for each bin pair.
			out <- .Call(cxx_count_background, pairs, new.pts$id, flank, exclude, filter, 
				bin.by.chr$first[[anchor2]], bin.by.chr$last[[anchor2]], bin.by.chr$first[[anchor1]], bin.by.chr$last[[anchor1]])
			if (is.character(out)) { stop(out) }
			if (!length(out[[1]])) { next }

			# Storing counts and locations. 
			if (any(out[[1]] < out[[2]])) { stop("anchor2 ID should not be less than anchor1 ID") }
			out.a[[idex]] <- out[[1]]
 			out.t[[idex]] <- out[[2]]
			out.counts[[idex]] <- out[[3]]

            # Adding the neighbourhood statistics.
            for (m in seq_along(modes)) {
                neighbor.counts[[m]][[idex]] <- out[[4]][[m]]
                neighbor.N[[m]][[idex]] <- out[[5]][[m]]
            }
			idex <- idex+1L
		}
	}

	# Collating all the other results.
	out.a <- unlist(out.a)
	out.t <- unlist(out.t)
	out.counts <- do.call(rbind, out.counts)

    all.assays <- list(counts=out.counts)
    for (m in seq_along(modes)) { 
        all.assays[[modes[m]]] <- do.call(rbind, neighbor.counts[[m]])
    }

	out.IS <- InteractionSet(all.assays, colData=DataFrame(totals=full.sizes), 
		interactions=GInteractions(anchor1=out.a, anchor2=out.t, regions=new.pts$region, mode="reverse"), 
        metadata=List(param=param, width=width))
   
    n.names <- .neighbor_numbers() 
    for (m in seq_along(modes)) { 
        mcols(out.IS)[[n.names[m]]] <- unlist(neighbor.N[[m]])
    }
    return(out.IS)
}
