getPairData <- function(file, param)
# This retrieves the fragment sizes, relative orientations and gaps from each directory produced by 
# preparePairs. This is a convenience function which allows people to avoid loading the entire directory 
# in (or manually parsing the said directory).
#
# written by Aaron Lun
# created 20 September 2014
# last modified 22 March 2017
{
    parsed <- .parseParam(param)
	allstuff <- .loadIndices(file, parsed$chrs)
	alll <- allo <- alli <- vector("list", sum(lengths(allstuff)))
	ix <- 1L

	# Running through all pairs.	
	for (ax in names(allstuff)) {
		current <- allstuff[[ax]] 
		for (tx in names(current)) { 
			extracted <- .getPairs(file, ax, tx)
			yielded <- .getStats(extracted, ax==tx, param$fragments)
			alll[[ix]] <- yielded$length
			allo[[ix]] <- yielded$orientation
			alli[[ix]] <- yielded$insert
			ix <- ix + 1L
		}
	}

	# Return objects.
	return(data.frame(length=unlist(alll), 
		orientation=unlist(allo), insert=unlist(alli)))
}

.getStats <- function(incoming, same.chr, fragments) 
# Gets the statistics for everything, including the fragment length,
# strand orientation and gap size.
{
	if (!nrow(incoming)) {
		return(list(length=integer(0), orientation=integer(0), gap=integer(0)))
	}
	output <- .Call(cxx_pair_stats, incoming$anchor1.id, incoming$anchor2.id, incoming$anchor1.pos, incoming$anchor2.pos,
		incoming$anchor1.len, incoming$anchor2.len, same.chr, start(fragments), end(fragments))
	if (is.character(output)) { stop(output) }
	names(output) <- c("length", "orientation", "insert")
	return(output)	
}

