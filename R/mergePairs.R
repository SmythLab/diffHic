mergePairs <- function(files, file.out)
# This function merges one or more separate count files together. It also produces a new index file
# representing the updated values of the older count file combination. The algorithm proceeds
# by reading all the index files in; splitting by anchor/target combinations, and then pulling out
# and combining the reads. This avoids having to read the entire structure into memory at once.
#
# written by Aaron Lun
# some time ago
# last modified 22 July 2015
{
	# Use a temporary file as a placeholder just in case 'file.out' is in 'files'.
	overall <- suppressWarnings(.loadIndices(files))
	tmpf <- tempfile(tmpdir=".")
	.initializeH5(tmpf) 
	on.exit({ if (file.exists(tmpf)) { unlink(tmpf, recursive=TRUE) } })

	# Merging for each combination, as necessary.
	for (ac in names(overall)) {
		current <- overall[[ac]]
		.addGroup(tmpf, ac)
		for (tc in names(current)) {
			fnames<-current[[tc]]

			out <- list()
			my.names <- NULL	
			for (ix in seq_along(fnames)) {
				if (!fnames[ix]) { next }
				xx <- .getPairs(files[ix], ac, tc)
				if (is.null(my.names)) {
					my.names <- colnames(xx) 
				} else if (!identical(my.names, colnames(xx))) {
					warning("column names are not identical between objects to be merged")
				}
				out[[length(out)+1L]] <- xx
			}

			# No need to protect against an empty list; there must be one non-empty element for .loadIndices to get here.
			out <- do.call(rbind, out)
			out <- out[order(out$anchor.id, out$target.id),]
			.writePairs(out, tmpf, ac, tc)
		}
	}

	# Moving the temporary, which is now the new file.
	if (!file.rename(tmpf, file.out)) { stop("cannot move file to the specified destination") }
	invisible(NULL)
}

