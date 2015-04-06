prunePairs <- function(file.in, param, file.out=file.in, max.frag=NA, min.inward=NA, min.outward=NA)
# This function collates all valid HiC pairs into a set of counts for each interaction. Basically, it 
# takes the directory output by 'prepareHiCPairs' and converts it into counts for each interaction. 
# If the fragment size is too big then we filter it out. We also filter out pairs which are too
# close together (min.*ward) with the difference between the two specifying inward/outward facing pairs.
#
# written by Aaron Lun
# created 9 September 2014
# last modified 20 March 2015
{
    # Use a temporary file as a placeholder, in case file.out==file.in.
	tmpf <- tempfile(tmpdir=".")
	on.exit({ if (file.exists(tmpf)) { unlink(tmpf) } })
	.initializeH5(tmpf)
	retained <- total <- by.len <- by.in <- by.out <- 0L

	# Parsing through the old index, counting/summing everything, and saving it to the
	# temporary file. We also remove any specified elements.
	allstuff <- .loadIndices(file.in, seqlevels(param$fragments))
	for (ax in names(allstuff)) { 
		current <- allstuff[[ax]]
		loaded <- FALSE
		for (tx in names(current)) { 
			collected <- .getPairs(file.in, ax, tx)
			stats <- .getStats(collected, ax==tx, param$fragments)

			if (!is.na(max.frag)) { 
				keep.len <- stats$length <= max.frag
				by.len <- by.len + sum(!keep.len)
			} else { keep.len <- TRUE }

			keep.in <- keep.out <- TRUE
			if (ax==tx) { 
				if (!is.na(min.inward)) { 
					keep.in <- stats$orientation!=1L | stats$insert >= min.inward
					by.in <- by.in + sum(!keep.in)
				}
				if (!is.na(min.outward)) { 
					keep.out <- stats$orientation!=2L | stats$insert >= min.outward
					by.out <- by.out + sum(!keep.out)
				}
			} 

			total <- total + nrow(collected) 
			collected <- collected[keep.len & keep.in & keep.out,]
			if (nrow(collected)) { 
				if (!loaded) { # Only adding a group if the data.frame is non-empty.
					.addGroup(tmpf, ax)
					loaded <- TRUE
				}
				.writePairs(collected, tmpf, ax, tx)
				retained <- retained + nrow(collected)
			}
		}
	}
	
	# Shuffling things around.
	if (!file.rename(tmpf, file.out)) { stop("cannot move file to the specified destination") }
	return(c(total=total, length=by.len, inward=by.in, outward=by.out, retained=retained))
}
