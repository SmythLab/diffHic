filterDirect <- function(data, ...)
# Implements the direct filtering method on the abundances of 
# inter-chromosomal bin pairs.
#
# written by Aaron Lun
# created 5 March 2015
# last modified 20 March 2015
{
	all.chrs <- seqnames(regions(data))
	is.inter <- as.logical(all.chrs[anchors(data, id=TRUE)]!=all.chrs[targets(data, id=TRUE)])
	ave.ab <- scaledAverage(asDGEList(data), ...)

	threshold <- .getInterThreshold(all.chrs, ave.ab[is.inter], empty=.makeEmpty(data, ...))
	return(list(abundances=ave.ab, threshold=threshold))
}

.getInterThreshold <- function(all.chrs, inter.ab, empty=NA) { 
	# Getting the total number of inter-chromosomal bins.
	n.bins <- as.numeric(runLength(all.chrs))
	total.bins <- sum(n.bins)
	n.inter <- total.bins * (total.bins + 1L)/2L - sum(n.bins * (n.bins + 1L)/2L)
	prop.kept <- length(inter.ab)/n.inter

	# Getting the threshold.
	if (prop.kept >= 1) { 
		threshold <- median(inter.ab) 
	} else if (prop.kept < 0.5) { 
		threshold <- empty
	} else { 
		threshold <- quantile(inter.ab, 1-0.5/prop.kept)
		names(threshold) <- NULL
	}
	return(threshold)
}

.makeEmpty <- function(data, ...) { scaledAverage(DGEList(0, lib.size=mean(data$totals)), ...) }

filterTrended <- function(data, span=0.25, ...)
# Implements the trended filtering method on the abundances of 
# inter-chromosomal bin pairs. 
#
# written by Aaron Lun
# created 5 March 2015
# last modified 20 March 2015
{
	dist <- getDistance(data, type="mid")
	log.dist <- log10(dist + exptData(data)$width)
	ave.ab <- scaledAverage(asDGEList(data), ...)

	# Filling in the missing parts of the interaction space.
	empty <- .makeEmpty(data, ...)
	is.intra <- !is.na(log.dist)
	n.intras <- sum(is.intra)
	all.chrs <- seqnames(regions(data))
	n.bins <- as.numeric(runLength(all.chrs))
	if (sum(n.bins * (n.bins + 1L)/2L) > 2*n.intras) { 
 	   	warning("too many missing regions in the intra-chromosomal interaction space to fill in") 
		trend.threshold <- loessFit(x=log.dist, y=ave.ab, span=span)$fitted
	} else {
		a.pts <- anchors(data, id=TRUE)[is.intra]
		t.pts <- targets(data, id=TRUE)[is.intra]
		o <- order(a.pts, t.pts) 
		a.pts <- a.pts[o]
		t.pts <- t.pts[o]

		extra.dist <- .Call(cxx_get_missing_dist, cumsum(runLength(all.chrs)),
			a.pts-1L, t.pts-1L, (start(regions(data))+end(regions(data)))/2)
		if (is.character(extra.dist)) { stop(extra.dist) }
		extra.dist <- log10(extra.dist + exptData(data)$width)
		trend.threshold <- loessFit(x=c(log.dist, extra.dist), 
			y=c(ave.ab, rep(empty, length(extra.dist))), 
			span=span)$fitted[1:length(log.dist)]
	}

	# Using the direct threshold.
	direct.threshold <- .getInterThreshold(seqnames(regions(data)), ave.ab[is.na(log.dist)], empty=empty)
	trend.threshold[is.na(log.dist)] <- direct.threshold
	return(list(abundances=ave.ab, threshold=trend.threshold, log.distance=log.dist)) 
}

