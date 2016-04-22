savePairs <- function(x, file, param)
# This function saves all counts in 'x' into a set of gzipped files in 'dir', along with
# an index file specifying the identity of each observed chromosome combination corresponding
# to each file. This speeds up any attempt at random access. The idea is to act as a 
# convenience function if you have a matrix of counts (or whatever) and you just want to save it.
#
# written by Aaron Lun 
# created some time ago
# last modified 22 April 2016
{
    if (any(colnames(x) %in% c("anchor.id", "target.id"))) { 
        stop("colnames 'anchor.*' and 'target.*' should be changed to 'anchor1.*' and 'anchor2.*'")
    }
	swap <- x$anchor1.id < x$anchor2.id
	if (any(swap)) { 
		temp <- x$anchor2.id[swap]
		x$anchor2.id[swap] <- x$anchor1.id[swap]
		x$anchor1.id[swap] <- temp
	}
	if (file.exists(file)) { unlink(file, recursive=TRUE) }

	# Need to reorder so fragments are sorted by chromosome COMBINATION.
	frag.out <- .splitByChr(param$fragments)
	all.chrs <- frag.out$chr
	full.chrs <- rep(seq_along(all.chrs), frag.out$last-frag.out$first+1L)
	achr <- full.chrs[x$anchor1.id]
	tchr <- full.chrs[x$anchor2.id]
	new.o <- order(achr, tchr, x$anchor1.id, x$anchor2.id)
	x <- x[new.o,]

	# Identifying stretches with the same chromatin pairs.
	new.achr <- achr[new.o]
	new.tchr <- tchr[new.o]
	if (length(new.achr) > 0L) {
		is.diff <- c(TRUE, diff(new.achr)!=0L | diff(new.tchr)!=0L)
		first.in.combo <- which(is.diff)
		last.in.combo <- c(first.in.combo[-1]-1L, length(new.o))
	} else { first.in.combo <- last.in.combo <- integer(0) }

	# Saving results.
	.initializeH5(file)
	for (ax in unique(new.achr[first.in.combo])) { .addGroup(file, all.chrs[ax]) }
	for (y in seq_along(first.in.combo)) {
		current <- first.in.combo[y]:last.in.combo[y]
		cur.a <- all.chrs[new.achr[current[1]]] 
		cur.t <- all.chrs[new.tchr[current[1]]]
	    .writePairs(x[current,], file, cur.a, cur.t)
	}
	invisible(NULL)
}
