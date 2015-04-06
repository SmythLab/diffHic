consolidatePairs <- function(data.list, result.list, equiweight=TRUE, combine.args=list()) 
# Consolidates results from multiple bin sizes. Returns a DIList of the
# relevant bin pairs in which things are now nested, as well as a table of
# combined p-values.
#
# written by Aaron Lun
# created 9 March 2015
{
	
	nset <- length(data.list)
	if (nset!=length(result.list)) { stop("data list must have same length as result list") }
	nall <- 0L
	for (x in 1:nset) {
		currows <- nrow(data.list[[x]])
		ntab <- nrow(result.list[[x]])
		if (currows!=ntab) { stop("corresponding entries of data and result lists must have same number of entries") }
		nall <- nall + ntab
	}

	# Boxing them and getting their weights.
	boxed <- do.call(boxPairs, data.list)
	if (equiweight) {
		weights <- list()
		for (x in 1:nset) { weights[[x]] <- 1/counts(boxed$pairs)[boxed$indices[[x]],x] }
		weights <- unlist(weights)
	} else {
		weights <- NULL
	}

	# Combining the p-values.
	all.tab <- do.call(rbind, result.list)
	result.com <- combineTests(ids=unlist(boxed$indices), tab=all.tab, weight=weights)
	return(list(id=boxed$indices, pairs=boxed$pairs, table=result.com))
}
