consolidatePairs <- function(indices, result.list, equiweight=TRUE, combine.args=list()) 
# Consolidates results from multiple bin sizes. Returns a DIList of the
# relevant bin pairs in which things are now nested, as well as a table of
# combined p-values.
#
# written by Aaron Lun
# created 9 March 2015
# last modified 22 July 2015
{
	
	nset <- length(indices)
	if (nset!=length(result.list)) { stop("indices must have same length as result list") }
	for (x in seq_len(nset)) {
		if (!identical(length(indices[[x]]), nrow(result.list[[x]]))) {
 		   	stop("corresponding entries of data and result lists must have same number of entries") }
	}

	# Boxing them and getting their weights.
	if (equiweight) {
		weights <- list()
		for (x in seq_len(nset)) {
			freq <- tabulate(indices[[x]])
			weights[[x]] <- 1/freq[indices[[x]]] 
		}
		weights <- unlist(weights)
	} else {
		weights <- NULL
	}

	# Combining the p-values.
	all.tab <- do.call(rbind, result.list)
	result.com <- do.call(combineTests, c(list(ids=unlist(indices), 
		tab=all.tab, weight=weights), combine.args))
	return(result.com)
}
