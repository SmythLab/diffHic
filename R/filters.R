filterDirect <- function(data, prior.count=2, reference=NULL)
# Implements the direct filtering method on the abundances of 
# inter-chromosomal bin pairs. Also allows for specification of
# a reference set of bin pairs (usually larger bins from which
# the abundances can be more stably computed).
#
# written by Aaron Lun
# created 5 March 2015
# last modified 16 November 2015
{
    .check_StrictGI(data)

	if (!is.null(reference)) { 
		stopifnot(identical(reference$totals, data$totals))
		scaling <- (.getBinSize(reference)/.getBinSize(data))^2
        ref <- .direct_filter(reference, prior.count=prior.count, scaling=scaling)
        actual.ab <- scaledAverage(asDGEList(data), prior.count=prior.count, scale=1)
		return(list(abundances=actual.ab, threshold=ref$threshold, ref=ref))
	}

    .direct_filter(data, prior.count=prior.count, scaling=1)
}

.direct_filter <- function(data, prior.count, scaling) 
# Actually does all the work of the direct filtering strategy.
# It computes the abundances and identifies the median for inter-chromosomal bin pairs.
{
   	all.chrs <- seqnames(regions(data))
	is.inter <- !intrachr(data)
	ave.ab <- scaledAverage(asDGEList(data), prior.count=prior.count, scale=scaling)
    empty.ab <- .makeEmpty(data, prior.count=prior.count, scale=scaling) 
	threshold <- .getInterThreshold(all.chrs, ave.ab[is.inter], empty=empty.ab)
    return(list(abundances=ave.ab, threshold=threshold))
}

.getBinSize <- function(data) 
# Gets the bin size in base pairs. This should be easy for bin pairs,
# but we also allow for more exotic set-ups, e.g., Capture-C loaded
# with regionCounts where anchors are bins around probes (evenly
# sized so treatable as bin pairs, but irregularly spaced).
{
	out <- metadata(data)$width
	if (is.null(out)) { out <- median(width(regions(data))) }
	return(out) 
}

.getInterThreshold <- function(all.chrs, inter.ab, empty=NA) 
# Computes the threshold from inter-chromosomal interactions.
# First we get the total number of inter-chromosomal bins,
# and then we compue the median (accounting for those lost).
{ 
	n.bins <- as.numeric(runLength(all.chrs))
	total.bins <- sum(n.bins)
	n.inter <- total.bins * (total.bins + 1L)/2L - sum(n.bins * (n.bins + 1L)/2L)
	n.inter <- max(n.inter, 1L) # avoid breakage when no inter-chromosomals are available.
	prop.kept <- length(inter.ab)/n.inter

	if (prop.kept >= 1) { 
		threshold <- median(inter.ab) 
	} else if (prop.kept < 0.5) { 
		warning("insufficient inter-chromosomal pairs for reliable threshold estimation")
		threshold <- empty
	} else { 
		threshold <- quantile(inter.ab, 1-0.5/prop.kept)
		names(threshold) <- NULL
	}
	return(threshold)
}

.makeEmpty <- function(data, ...) { 
    if (nrow(data)) { 
        y <- asDGEList(data[1,])
        y$counts[] <- 0L
    } else {
        y <- asDGEList(data)
        y$counts <- rbind(integer(ncol(data)))
    }
    scaledAverage(y, ...) 
}

filterTrended <- function(data, span=0.25, prior.count=2, reference=NULL)
# Implements the trended filtering method on the abundances of 
# inter-chromosomal bin pairs. Again, with allowances for a reference set.
#
# written by Aaron Lun
# created 5 March 2015
# last modified 16 November 2016
{
    .check_StrictGI(data)

	if (!is.null(reference)) {
        stopifnot(identical(reference$totals, data$totals))
        scaling <- (.getBinSize(reference)/.getBinSize(data))^2
        ref <- .trended_filter(reference, span=span, prior.count=prior.count, scaling=scaling)

        actual.ab <- scaledAverage(asDGEList(data), prior.count=prior.count, scale=1)
		actual.dist <- log10(pairdist(data, type="mid") + .getBinSize(data))
		
		new.threshold <- approx(x=ref$log.distance, y=ref$threshold, xout=actual.dist, rule=2)$y
		new.threshold[is.na(actual.dist)] <- ref$threshold[is.na(ref$log.distance)][1] # Direct threshold.

        return(list(abundances=actual.ab, threshold=new.threshold, log.distance=actual.dist, ref=ref)) 
	}

    .trended_filter(data, span=span, prior.count=prior.count, scaling=1)        
}

.trended_filter <- function(data, span, prior.count, scaling) {
	dist <- pairdist(data, type="mid")
	log.dist <- log10(dist + .getBinSize(data)) # Adding an appropriate prior.
	ave.ab <- scaledAverage(asDGEList(data), prior.count=prior.count, scale=scaling)

	# Filling in the missing parts of the interaction space.
	empty <- .makeEmpty(data, prior.count=prior.count, scale=scaling)
	is.intra <- !is.na(log.dist)
	n.intras <- sum(is.intra)
	all.chrs <- seqnames(regions(data))
	n.bins <- as.numeric(runLength(all.chrs))

    if (sum(n.bins * (n.bins + 1L)/2L) > 2*n.intras) { 
 	   	warning("too many missing regions in the intra-chromosomal interaction space to fill in") 
		trend.threshold <- loessFit(x=log.dist, y=ave.ab, span=span)$fitted
	} else {
		a.pts <- anchors(data, type="first", id=TRUE)[is.intra]
		t.pts <- anchors(data, type="second", id=TRUE)[is.intra]
		o <- order(a.pts, t.pts) 
		a.pts <- a.pts[o]
		t.pts <- t.pts[o]

		extra.dist <- .Call(cxx_get_missing_dist, cumsum(runLength(all.chrs)),
			a.pts-1L, t.pts-1L, (start(regions(data))+end(regions(data)))/2)
		if (is.character(extra.dist)) { stop(extra.dist) }
		extra.dist <- log10(extra.dist + .getBinSize(data))
		trend.threshold <- loessFit(x=c(log.dist, extra.dist), 
			y=c(ave.ab, rep(empty, length(extra.dist))), 
			span=span)$fitted[seq_along(log.dist)]
	}

	# Using the direct threshold.
	is.inter <- is.na(dist)
	if (any(is.inter)) { 
		direct.threshold <- .getInterThreshold(seqnames(regions(data)), ave.ab[is.inter], empty=empty)
		trend.threshold[is.inter] <- direct.threshold
	}
	return(list(abundances=ave.ab, threshold=trend.threshold, log.distance=log.dist)) 
}

filterDiag <- function(data, by.dist=0, by.diag=0L, dist, ...)
# Filters diagonal elements, with options for supplying
# your own distance if you've already computed it.
# Can be extended to near-diagonal elements with by.dist or by.diag,
# which compute based on distance and bin-based diagonal, respectively.
#
# written by Aaron Lun
# created 6 October 2015
# last modified 8 December 2015
{
    .check_StrictGI(data)

	if (missing(dist)) { 
		dist <- pairdist(data, ...) 
	} else {
		dist <- as.numeric(dist)
		stopifnot(length(dist)==nrow(data)) 
	}
	by.dist <- as.numeric(by.dist)
	keep <- dist > by.dist 

	by.diag <- as.integer(by.diag)
	if (by.diag) { 
		diag.level <- pairdist(data, type="diag")
		keep <- keep & diag.level > by.diag
	}

	keep <- keep | is.na(dist) 
	return(keep)
}
