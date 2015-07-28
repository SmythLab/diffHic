correctedContact <- function(data, iterations=50, exclude.local=1, ignore.low=0.02, winsor.high=0.02, 
	average=TRUE, dist.correct=FALSE)
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
# last modified 22 July 2015
{
	if (!average & ncol(data)>1L) {
		collected.truth <- collected.bias <- collected.max <- collected.trend <- list()
		for (lib in seq_len(ncol(data))) {
			out <- Recall(data[,lib], iterations=iterations, exclude.local=exclude.local, ignore.low=ignore.low, 
				winsor.high=winsor.high, average=FALSE, dist.correct=dist.correct)
			collected.truth[[lib]] <- out$truth
			collected.bias[[lib]] <- out$bias
			collected.max[[lib]] <- out$max
			collected.trend[[lib]] <- out$trend
		}

		output <- list(truth=do.call(cbind, collected.truth), 
			bias=do.call(cbind, collected.bias),
			max=do.call(cbind, collected.max))
		if (dist.correct) { output$trend <- do.call(cbind, collected.trend) }
		return(output) 
	}

	# Checking arguments.
	iterations <- as.integer(iterations)
	if (iterations <= 0L) { stop("number of iterations must be a positive integer") }
	ignore.low <- as.double(ignore.low)
	if (ignore.low >= 1) { stop("proportion of low coverage fragments to ignore should be less than 1") }
  	winsor.high <- as.double(winsor.high)
	if (winsor.high >= 1) { stop("proportion of high coverage interactions to winsorize should be less than 1") }
	exclude.local <- as.integer(exclude.local)
    
	# Computing average counts, with or without distance correction.
	if (dist.correct) { 
		temp <- data
		temp$totals <- 1e6 # need library size to be reflected in fitted value of trend.
		trended <- filterTrended(temp, prior.count=0)
		ave.counts <- 2^(trended$abundance - trended$threshold)
		is.local <- !is.na(trended$log.distance)
		nzero <- !is.na(ave.counts)
	} else {
		is.local <- !is.na(getDistance(data))
   		log.lib <- log(data$totals)
		if (length(log.lib)>1L) {
			ave.counts <- exp(edgeR::mglmOneGroup(counts(data), offset=log.lib - mean(log.lib)))
			nzero <- !is.na(ave.counts)
		} else {
			nzero <- counts(data) != 0L
			ave.counts <- as.double(counts(data))
		}
	}

	out<-.Call(cxx_iterative_correction, ave.counts[nzero], anchors(data, id=TRUE)[nzero], targets(data, id=TRUE)[nzero], 
		is.local[nzero], length(regions(data)), iterations, exclude.local, ignore.low, winsor.high)
 	if (is.character(out)) { stop(out) }
	full.truth <- rep(0, length(nzero))
	full.truth[nzero] <- out[[1]]
	out[[1]] <- full.truth
	
	names(out) <- c("truth", "bias", "max")
	if (dist.correct) { out$trend <- trended$threshold }
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

