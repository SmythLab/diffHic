clusterPairs <- function(..., tol, upper=1e6) 
# This function examines the bin pairs in two-dimensional space and 
# clusters those pairs which are close together. Specifically, it does 
# so if the Chebyshev distance between the regions is less than 'tol'.
# It also splits clusters which are greater than 'upper' by partitioning
# them into smaller clusters of equal size.
#
# written by Aaron Lun
# created 6 December 2013
# last modified 4 June 2015
{
	tol <- as.integer(tol)
	stopifnot(tol>=0L) # Minimum overlap not supported.
	upper <- as.integer(upper)

	all.data <- list(...)
	achrs <- tchrs <- astarts <- aends <- tstarts <- tends <- list()
	for (x in 1:length(all.data)) { 
		data <- all.data[[x]]
		region <- regions(data)
		allchrs <- as.character(seqnames(region))
		aid <- anchors(data, id=TRUE)
		tid <- targets(data, id=TRUE)

		achrs[[x]] <- allchrs[aid]
		tchrs[[x]] <- allchrs[tid]
		astarts[[x]] <- start(region)[aid]
		aends[[x]] <- end(region)[aid]+1L
		tstarts[[x]] <- start(region)[tid]
		tends[[x]] <- end(region)[tid]+1L
	}

	achrs <- unlist(achrs)
	tchrs <- unlist(tchrs)
	astarts <- unlist(astarts)
	tstarts <- unlist(tstarts)
	aends <- unlist(aends)
	tends <- unlist(tends)

	# Figuring out which are blocks of separate chromosomes.	
	ro <- order(achrs, tchrs)
	achrs <- achrs[ro]
	tchrs <- tchrs[ro]
	astarts <- astarts[ro]
	tstarts <- tstarts[ro]
	aends <- aends[ro]
	tends <- tends[ro]
	n <- length(ro)
	is.new <- which(c(TRUE, achrs[-1]!=achrs[-n] | tchrs[-1]!=tchrs[-n]))
	upnext <- c(is.new[-1]-1L, n)

	# Now, running through.
	all.ids <- integer(n)
	bonus <- 0L
	for (i in 1:length(upnext)) {
		current <- is.new[i]:upnext[i]			
		curas <- astarts[current]	
		curae <- aends[current]	
		curts <- tstarts[current]	
		curte <- tends[current]
		po <- order(curas, curts)
	
		out <- .Call(cxx_cluster_2d, curas[po], curts[po], curae[po], curte[po], tol, FALSE)
		if (is.character(out)) { stop(out) }
		out[po] <- out
		if (length(upper)) { 
			out <- .Call(cxx_split_clusters, out, curas, curts, curae, curte, upper)
			if (is.character(out)) { stop(out) }
			xo <- order(out)
			out[xo]<-cumsum(c(TRUE, diff(out[xo])!=0L)) # Cleaning it up a little, to keep the IDs reasonably tight.
		}

		all.ids[current] <- out + bonus
		bonus <- bonus + max(out)
	}

	# Getting the bounding box for each cluster.
	min.box <- .minBoundingBox(all.ids, achrs, astarts, aends, tchrs, tstarts, tends, seqinfo(region))
	all.ids[ro] <- all.ids	
	indices <- list()
	last <- 0
	for (x in 1:length(all.data)) { 
		currows <- nrow(all.data[[x]])
		if (!currows) { next } 
		indices[[x]] <- all.ids[last + 1:currows]
		last <- last + currows
	}
	names(indices) <- names(all.data)
	return(list(indices=indices, anchors=min.box$anchors, targets=min.box$targets))
}

# No need to consider special behaviour beyond the diagonal.
# If a bin pair is within Chebyshev distance of a target bin pair but beyond the diagonal,
# then its reflection will also be within that distance to the target bin pair.
# It's like trying to reflect a corner of a square around the diagonal, the reflection will just lie in the square itself.

.minBoundingBox <- function(all.ids, achrs, astarts, aends, tchrs, tstarts, tends, seqinf) {
	a.out <- .Call(cxx_get_bounding_box, all.ids, astarts, aends)
	if (is.character(a.out)) { stop(a.out) }
	anchor.bounds <- GRanges(achrs[a.out[[1]]], IRanges(a.out[[2]], a.out[[3]] - 1L), seqinfo=seqinf)
	t.out <- .Call(cxx_get_bounding_box, all.ids, tstarts, tends)
	if (is.character(t.out)) { stop(t.out) }
	target.bounds <- GRanges(tchrs[t.out[[1]]], IRanges(t.out[[2]], t.out[[3]] - 1L), seqinfo=seqinf)
	return(list(anchors=anchor.bounds, targets=target.bounds))
}
