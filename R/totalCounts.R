totalCounts <- function(files, param)
# This function gets the total counts in a bunch of files.  This is designed
# for whenever the total counts must be rapidly extracted, without the need to
# count across the interaction space.
#
# written by Aaron Lun
# created 17 September 2014
# last modified 20 March 2015
{
	nlibs <- length(files)
	if (nlibs==0) { 
		stop("number of libraries must be positive")
	} 
	fragments <- param$fragments
 	chrs <- seqlevels(fragments) 
	full.sizes <- integer(nlibs)

	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
    for (anchor in names(overall)) {
        current <- overall[[anchor]]
		for (target in names(current)) {
			if (!.checkIfPairOK(restrict, anchor, target)) { next }

			# Getting totals.
			pairs <- .baseHiCParser(current[[target]], files, anchor, target, discard=discard, cap=cap)
			for (lib in 1:length(pairs)) { full.sizes[lib] <- full.sizes[lib] + nrow(pairs[[lib]]) }
		}
	}

	return(full.sizes)
}


