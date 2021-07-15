totalCounts <- function(files, param)
# This function gets the total counts in a bunch of files.  This is designed
# for whenever the total counts must be rapidly extracted, without the need to
# count across the interaction space.
#
# written by Aaron Lun
# created 17 September 2014
# last modified 14 May 2017
{
	nlibs <- length(files)
	if (nlibs==0L) { stop("number of libraries must be positive") }
	full.sizes <- integer(nlibs)

	# Running through each pair of chromosomes.
	loadfuns <- preloader(files, param=param)
    for (anchor in names(loadfuns)) {
        current <- loadfuns[[anchor]]
		for (target in names(current)) {
            curfuns <- current[[target]]

			# Getting totals.
            for (lib in seq_len(nlibs)) { 
                libfun <- curfuns[[lib]]
                if (!is.null(libfun)) {
                    full.sizes[lib] <- .addToTotal(full.sizes[lib], nrow(libfun()))
                }
            }
		}
	}

	return(full.sizes)
}


