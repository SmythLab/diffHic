filterDirect <- function(data, prior.count=2, reference=NULL)
# Implements the direct filtering method on the abundances of 
# inter-chromosomal bin pairs. Also allows for specification of
# a reference set of bin pairs (usually larger bins from which
# the abundances can be more stably computed).
#
# written by Aaron Lun
# created 5 March 2015
# last modified 24 June 2015
{
	if (!is.null(reference)) { 
		actual.ab <- scaledAverage(asDGEList(data), prior.count=prior.count)
		ref <- Recall(reference, prior.count=prior.count)

		stopifnot(identical(reference$totals, data$totals))
		scaling <- (.getBinSize(reference)/.getBinSize(data))^2
		adj.thresh <- .repriorAveLogCPM(ref$threshold, totals=data$totals,
			prior.count=prior.count, scaling=scaling)
		return(list(abundances=actual.ab, threshold=adj.thresh, ref=ref))
	}

	all.chrs <- seqnames(regions(data))
	is.inter <- as.logical(all.chrs[anchors(data, id=TRUE)]!=all.chrs[targets(data, id=TRUE)])
	ave.ab <- scaledAverage(asDGEList(data), prior.count=prior.count)

	threshold <- .getInterThreshold(all.chrs, ave.ab[is.inter],
		empty=.makeEmpty(data, prior.count=prior.count))
	return(list(abundances=ave.ab, threshold=threshold))
}

.getBinSize <- function(data) 
# Gets the bin size in base pairs. This should be easy for bin pairs,
# but we also allow for more exotic set-ups, e.g., Capture-C loaded
# with regionCounts where anchors are bins around probes (evenly
# sized so treatable as bin pairs, but irregularly spaced).
{
	out <- exptData(data)$width
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

.makeEmpty <- function(data, ...) { scaledAverage(DGEList(rbind(integer(ncol(data))), lib.size=data$totals), ...) }

.repriorAveLogCPM <- function(AveLogCPM, totals, prior.count, scaling)
# Adjusting the average log-CPM to use a new prior count.
{
	ave.count <- 2^AveLogCPM * mean(totals) / 1e6
	ave.count <- ave.count + prior.count*(scaling - 1)
	return(log2(ave.count * 1e6 / mean(totals) / scaling))
}

filterTrended <- function(data, span=0.25, prior.count=2, reference=NULL)
# Implements the trended filtering method on the abundances of 
# inter-chromosomal bin pairs. Again, with allowances for a reference set.
#
# written by Aaron Lun
# created 5 March 2015
# last modified 22 July 2015
{
	if (!is.null(reference)) {
		actual.ab <- scaledAverage(asDGEList(data), prior.count=prior.count)
		actual.dist <- log10(getDistance(data, type="mid") + .getBinSize(data))
		ref <- Recall(reference, span=span, prior.count=prior.count)
		
		new.threshold <- approx(x=ref$log.distance, y=ref$threshold, xout=actual.dist, rule=2)$y
		new.threshold[is.na(actual.dist)] <- ref$threshold[is.na(ref$log.distance)][1] # Direct threshold.

		stopifnot(identical(reference$totals, data$totals))
		scaling <- (.getBinSize(reference)/.getBinSize(data))^2
		adj.thresh <- .repriorAveLogCPM(new.threshold, totals=data$totals,
			prior.count=prior.count, scaling=scaling)
		return(list(abundances=actual.ab, threshold=adj.thresh, log.distance=actual.dist, ref=ref)) 
	}

	dist <- getDistance(data, type="mid")
	log.dist <- log10(dist + .getBinSize(data))
	ave.ab <- scaledAverage(asDGEList(data), prior.count=prior.count)

	# Filling in the missing parts of the interaction space.
	empty <- .makeEmpty(data, prior.count=prior.count)
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
		extra.dist <- log10(extra.dist + .getBinSize(data))
		trend.threshold <- loessFit(x=c(log.dist, extra.dist), 
			y=c(ave.ab, rep(empty, length(extra.dist))), 
			span=span)$fitted[seq_along(log.dist)]
	}

	# Using the direct threshold.
	is.inter <- is.na(log.dist)
	if (any(is.inter)) { 
		direct.threshold <- .getInterThreshold(seqnames(regions(data)), ave.ab[is.inter], empty=empty)
		trend.threshold[is.inter] <- direct.threshold
	}
	return(list(abundances=ave.ab, threshold=trend.threshold, log.distance=log.dist)) 
}

