clusterPairs <- function(data, tol, upper=1e6) 
# This function examines the bin pairs in two-dimensional space and 
# clusters those pairs which are close together. Specifically, it does 
# so if the Chebyshev distance between the regions is less than 'tol'.
# It also splits clusters which are greater than 'upper' by partitioning
# them into smaller clusters of equal size.
#
# written by Aaron Lun
# created 6 December 2013
# last modified 24 March 2015
{
	region <- regions(data)
	allchrs <- as.character(seqnames(region))
	aid <- anchors(data, id=TRUE)
	tid <- targets(data, id=TRUE)

	achrs <- allchrs[aid]
	tchrs <- allchrs[tid]
	tol <- as.integer(tol)
	stopifnot(tol>=0L) # Minimum overlap not supported.
	upper <- as.integer(upper)

	# Figuring out which are blocks of separate chromosomes.	
	ro <- order(achrs, tchrs)
	achrs <- achrs[ro]
	tchrs <- tchrs[ro]
	n <- length(ro)
	is.new <- which(c(TRUE, achrs[-1]!=achrs[-n] | tchrs[-1]!=tchrs[-n]))
	upnext <- c(is.new[-1]-1L, n)
	
	# Setting up the starts and ends.
	astarts <- start(region[aid[ro]])
	aends <- end(region[aid[ro]])+1L
	tstarts <- start(region[tid[ro]])
	tends <- end(region[tid[ro]])+1L

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
	
		out <- .Call(cxx_cluster_2d, curas[po], curts[po], curae[po], curte[po], tol)
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
	a.out <- .Call(cxx_get_bounding_box, all.ids, astarts, aends)
	if (is.character(a.out)) { stop(a.out) }
	anchor.bounds <- GRanges(achrs[a.out[[1]]], IRanges(a.out[[2]], a.out[[3]] - 1L), seqinfo=seqinfo(region))
	t.out <- .Call(cxx_get_bounding_box, all.ids, tstarts, tends)
	if (is.character(t.out)) { stop(t.out) }
	target.bounds <- GRanges(tchrs[t.out[[1]]], IRanges(t.out[[2]], t.out[[3]] - 1L), seqinfo=seqinfo(region))
	
	all.ids[ro] <- all.ids
	return(list(id=all.ids, anchors=anchor.bounds, targets=target.bounds))
}

# No need to consider special behaviour beyond the diagonal.
# If a bin pair is within Chebyshev distance of a target bin pair but beyond the diagonal,
# then its reflection will also be within that distance to the target bin pair.
# It's like trying to reflect a corner of a square around the diagonal, the reflection will just lie in the square itself.

