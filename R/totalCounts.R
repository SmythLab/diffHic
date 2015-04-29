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
	frag.by.chr <- .splitByChr(fragments)
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

			# Getting totals.
			pairs <- .baseHiCParser(current[[target]], files, anchor, target, 
				chr.limits=frag.by.chr, discard=discard, cap=cap)
			full.sizes <- full.sizes + sapply(pairs, FUN=nrow)
		}
	}

	return(full.sizes)
}


