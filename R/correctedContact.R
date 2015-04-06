correctedContact <- function(data, iterations=50, exclude.local=1, ignore.low=0.02, winsor.high=0.02, average=TRUE, dispersion=0.05)
# This performs the iterative correction method of Mirny et al. (2012) to
# identify the true contact probability of each patch of the interaction
# space. The idea is to use the true contact probability as a filter
# to identify the top set of patches (assuming that most patches represent
# non-specific ligation). Note that for this to work, you need 'filter=1L' 
# during count loading. There is no alternative to loading everything in,
# because you need to compute the truth by dividing and resumming each count.
#
# written by Aaron Lun
# some time ago	
# last modified 20 March 2015
{
	if (!average & ncol(data)>1L) {
		collected.truth <- collected.bias <- collected.max <- list()
		for (lib in 1:ncol(data)) { 
			out <- Recall(data[,lib], iterations=iterations, exclude.local=exclude.local, ignore.low=ignore.low, 
				winsor.high=winsor.high, average=FALSE, dispersion=dispersion)
			collected.truth[[lib]] <- out$truth
			collected.bias[[lib]] <- out$bias
			collected.max[[lib]] <- out$max
		}
		return(list(truth=do.call(cbind, collected.truth), 
			bias=do.call(cbind, collected.bias),
			max=do.call(cbind, collected.max)))
	}

	# Checking arguments.
	iterations <- as.integer(iterations)
	if (iterations <= 0L) { stop("number of iterations must be a positive integer") }
	ignore.low <- as.double(ignore.low)
	if (ignore.low >= 1) { stop("proportion of low coverage fragments to ignore should be less than 1") }
  	winsor.high <- as.double(winsor.high)
	if (winsor.high >= 1) { stop("proportion of high coverage interactions to winsorize should be less than 1") }
	exclude.local <- as.integer(exclude.local)
    
	is.local <- !is.na(getDistance(data))
   	log.lib <- log(data$totals)
	if (length(log.lib)>1L) {
		ave.counts <- exp(edgeR::mglmOneGroup(counts(data), offset=log.lib - mean(log.lib), dispersion=dispersion))
		nzero <- !is.na(ave.counts)
	} else {
		nzero <- counts(data) != 0L
		ave.counts <- as.double(counts(data))
	}

	out<-.Call(cxx_iterative_correction, ave.counts[nzero], anchors(data, id=TRUE)[nzero], targets(data, id=TRUE)[nzero], 
		is.local[nzero], length(regions(data)), iterations, exclude.local, ignore.low, winsor.high)
 	if (is.character(out)) { stop(out) }
	full.truth <- rep(NA, length(nzero))
	full.truth[nzero] <- out[[1]]
	out[[1]] <- full.truth
	
	names(out) <- c("truth", "bias", "max")
	return(out)
}

# getBias <- function(bias, index, na.action=min)
# # As its name suggests, gets the bias. Some protection is provided against
# # NA values for the bias. These generally correspond to low-abundance regions 
# # for which iterative correction does not provide stable results. We replace 
# # it with the next-best thing, usually the minimum of the non-NA values.
# {
# 	output <- bias[index]
# 	if (any(is.na(output))) { output[is.na(output)] <- na.action(bias[!is.na(bias)]) }
# 	return(output)
# }

