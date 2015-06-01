compartmentalize <- function(data, centers=2, dist.correct=TRUE,
		cov.correct=TRUE, robust.cov=5, inter=FALSE, ...)
# Computes compartments for every chromosome, using the intra-chromosomal
# contact maps that have been corrected for distance effects.
#
# written by Aaron Lun
# created 26 May 2015
# last modified 28 May 2015
{
	if (!inter) {
		is.intra <- !is.na(getDistance(data))
		data <- data[is.intra,]
	}
	if (dist.correct) { 
		trended <- filterTrended(data)
		contacts <- trended$abundance - trended$threshold
		dist2trend <- approxfun(x=trended$log.distance, y=trended$threshold, rule=2)
	} else {
		contacts <- aveLogCPM(asDGEList(data))
		dist2trend <- function(x) { 0 } # Trend correction function does nothing if no distance correction is requested.
	}

	if (!inter) { 
		# Going chromosome-by-chromosome.
		stored <- list()
		for (chr in seqlevels(regions(data))) { 
			mat <- as.matrix(data, first=chr, fill=contacts)
			stored[[chr]] <- .compartChr(mat, data, dist2trend, robust.cov, cov.correct, centers, ...)
		}
	} else {
		# Using the entirety of the interaction space.				
		mat <- as.matrix(data, fill=contacts)
		stored <- .compartChr(mat, data, dist2trend, robust.cov, cov.correct, centers, ...)
	}

	return(stored)
}

.compartChr <- function(mat, data, dist2trend, robust.cov, cov.correct, centers, ...) {
	# Filling NA's (i.e., zero's). Using mid-distance, interpolating to get the trend.
	lost <- which(is.na(mat), arr.ind=TRUE)
	lost.dist <- abs(mid(ranges(anchors(data[lost[,1]]))) - mid(ranges(targets(data[lost[,2]]))))
	lost.dist <- log10(lost.dist + exptData(data)$width)
	mat[is.na(mat)] <- .makeEmpty(data) - dist2trend(lost.dist)

	# Correcting for coverage biases, by subtracting half the average coverage from both rows 
	# and columns. This is equivalent to dividing by square root of coverage, which works pretty 
	# well in place of a more rigorous iterative approach (check out Rao's supplementaries).
	rwm <- log2(rowMeans(2^mat))
	if (cov.correct) { 
		mat <- mat - rwm/2
		mat <- t(t(mat) - rwm/2)
	}
	
	# Robustifying.
	if (!is.na(robust.cov)) {
 	    rwm.med <- median(rwm)
		rwm.mad <- mad(rwm, center=rwm.med)	
		keep <- (rwm <= rwm.med + robust.cov*rwm.mad) & (rwm >= rwm.med - robust.cov*rwm.mad)
		temp.mat <- mat
		mat <- mat[keep,keep,drop=FALSE] 
	}

	# K-means.
	if (nrow(mat) > centers) {
		out <- kmeans(mat, centers=centers, ...)
		comp <- out$cluster
	} else {
		comp <- 1:nrow(mat)
	}

	# Filling the robustified bins back in.
	if (!is.na(robust.cov)) { 
		mat <- temp.mat
		temp <- integer(length(keep))
		temp[keep] <- comp
		temp[!keep] <- NA
		comp <- temp
	}

	names(comp) <- rownames(mat)
	return(list(compartment=comp, matrix=mat))
}


