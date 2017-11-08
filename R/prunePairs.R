prunePairs <- function(file.in, param, file.out=file.in, max.frag=NA, min.inward=NA, min.outward=NA)
# This function collates all valid HiC pairs into a set of counts for each interaction. Basically, it 
# takes the directory output by 'prepareHiCPairs' and converts it into counts for each interaction. 
# If the fragment size is too big then we filter it out. We also filter out pairs which are too
# close together (min.*ward) with the difference between the two specifying inward/outward facing pairs.
#
# written by Aaron Lun
# created 9 September 2014
# last modified 8 November 2017
{
    # Use a temporary file as a placeholder, in case file.out==file.in.
	tmpf <- tempfile(tmpdir=".")
	on.exit({ if (file.exists(tmpf)) { unlink(tmpf) } })
	.initializeH5(tmpf)
	retained <- total <- by.len <- by.in <- by.out <- 0L
   
    # Figuring out which chromosomes are of interest. 
    parsed <- .parseParam(param, bin=FALSE)
    chrs <- parsed$chrs
    frag.by.chr <- parsed$frag.by.chr
    restrict <- parsed$restrict
	allstuff <- .loadIndices(file.in, chrs, restrict)

    # Figuring out if filtering on the read pairs is necessary.
    discard <- parsed$discard
    cap <- parsed$cap
    if (.isDNaseC(param)) { 
        cap <- NA # no binning here, so cap is useless.
    }
   
	# Parsing through the old index.
	for (ax in names(allstuff)) { 
		current <- allstuff[[ax]]
		loaded <- FALSE
		for (tx in names(current)) { 
            collected <- .baseHiCParser(TRUE, file.in, ax, tx, chr.limits=frag.by.chr, 
                discard=discard, cap=NA, width=NA, retain=NULL)[[1]]
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
