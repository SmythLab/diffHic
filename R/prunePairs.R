prunePairs <- function(file.in, param, file.out=file.in, max.frag=NA, min.inward=NA, min.outward=NA)
# This function collates all valid HiC pairs into a set of counts for each interaction. Basically, it 
# takes the directory output by 'prepareHiCPairs' and converts it into counts for each interaction. 
# If the fragment size is too big then we filter it out. We also filter out pairs which are too
# close together (min.*ward) with the difference between the two specifying inward/outward facing pairs.
#
# written by Aaron Lun
{
    # Use a temporary file as a placeholder, in case file.out==file.in.
    # We use the same directory to avoid cross-device links.
	tmpf <- tempfile(tmpdir=dirname(file.out))
	on.exit({ if (file.exists(tmpf)) { unlink(tmpf) } })
	.initializeH5(tmpf)
	retained <- total <- by.len <- by.in <- by.out <- 0L

    # Capping is handled after extraction, and needs to be turned off in 'preloader'.
    cap <- param$cap
    if (.isDNaseC(param)) { 
        cap <- NA_integer_ # no restriction fragments or binning, so capping is useless.
    }
    param <- reform(param, cap=NA_integer_)

	# Parsing through the old index.
    loadfuns <- preloader(file.in, param=param, retain=NULL)
	for (ax in names(loadfuns)) { 
		current <- loadfuns[[ax]]
		loaded <- FALSE

		for (tx in names(current)) { 
            curfun <- current[[tx]][[1]]
            collected <- curfun() # can't be NULL, only one library!
			stats <- .getStats(collected, ax==tx, param$fragments)

            # Removing reads above the max.
			if (!is.na(max.frag)) { 
				keep.len <- stats$length <= max.frag
				by.len <- by.len + sum(!keep.len)
			} else { 
                keep.len <- TRUE 
            }

            # Removing inward/outward facing reads on the same chromosome.
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
            
            # Summing and filtering.
			total <- total + nrow(collected) 
            collected <- collected[keep.len & keep.in & keep.out,]

            # Removing reads above the cap for a given restriction fragment pair.
            # This is done AFTER other pruning for consistency with squareCounts(), had
            # we just performed the other pruning and used the cap in squareCounts().
            # The total needs to be adjusted as these reads are theoretically lost prior to 
            # the other pruning; this ensures that the total-retained is only due to the other pruning.
            if (!is.na(cap)) {
                capped <- .Call(cxx_cap_input, collected$anchor1.id, collected$anchor2.id, cap)
                collected <- collected[capped,]
                total <- total - sum(!capped) 
            }

            # Saving to a new file.
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
	if (!file.rename(tmpf, file.out)) { 
        stop("cannot move file to the specified destination") 
    }
	return(c(total=total, length=by.len, inward=by.in, outward=by.out, retained=retained))
}
