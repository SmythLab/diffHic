mergePairs <- function(files, file.out)
# This function merges one or more separate count files together. It also produces a new index file
# representing the updated values of the older count file combination. The algorithm proceeds
# by reading all the index files in; splitting by anchor combinations, and then pulling out
# and combining the reads. This avoids having to read the entire structure into memory at once.
#
# written by Aaron Lun
{
	# Use a temporary file (in the same directory as the output file) 
    # as a placeholder just in case 'file.out' is in 'files'.
	tmpf <- tempfile(tmpdir=dirname(file.out))
	.initializeH5(tmpf) 
	on.exit({ if (file.exists(tmpf)) { unlink(tmpf, recursive=TRUE) } })

	# Merging for each combination, as necessary.
    loadfuns <- preloader(files, param=NULL, retain=NULL)
	for (ac in names(loadfuns)) {
		current <- loadfuns[[ac]]
		.addGroup(tmpf, ac)
		for (tc in names(current)) {
			curfuns <- current[[tc]]

			out <- vector("list", length(curfuns)) 
			for (ix in seq_along(curfuns)) {
                curFUN <- curfuns[[ix]]
                if (is.null(curFUN)) { 
                    next 
                }
				out[[ix]] <- curFUN()
			}

			# No need to protect against an empty list; 
            # there must be one non-empty element for preloader to get here.
			out <- do.call(rbind, out)
			out <- out[order(out$anchor1.id, out$anchor2.id),]
			.writePairs(out, tmpf, ac, tc)
		}
	}

	# Moving the temporary, which is now the new file.
	if (!file.rename(tmpf, file.out)) { 
        stop("cannot move file to the specified destination") 
    }
	invisible(NULL)
}

