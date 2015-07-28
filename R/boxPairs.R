boxPairs <- function(..., reference, minbox=FALSE)
# This function reports bin pairs that are nested within other bin pairs.  The
# idea is to consolidate smaller bin pairs into their larger counterparts for
# summarization of analyses involving multiple bin sizes.
#
# written by Aaron Lun
# created 3 June 2014
# last modified 22 July 2015
{
	all.hits <- list(...)
	nk <- length(all.hits)
	if (missing(reference)) { 
		reference <- max(sapply(all.hits, FUN=function(x) { exptData(x)$width }))
	}
	fragments <- exptData(all.hits[[1]])$param$fragments
	for (x in all.hits[-1]) { 
		if (!identical(exptData(x)$param$fragments, fragments)) {
			stop("fragment boundaries should be the same between DIList objects")
		}
	}
	parents <- .getBinID(fragments, reference)$region

	# Collating all results in terms of parents.
	all.a <- all.t <- all.mode <- all.idx <- list()
	num.pairs <- list()
	seq.it.nk <- seq_len(nk)
	for (x in seq.it.nk) {
		current <- all.hits[[x]]
		ncur <- nrow(current)
		olap <- findOverlaps(regions(current), parents, type="within", select="first")
		if (any(is.na(olap))) { stop("smaller bins must be fully contained within larger bins") }
		
		all.a[[x]] <- olap[anchors(current, id=TRUE)]
		all.t[[x]] <- olap[targets(current, id=TRUE)]
		all.mode[[x]] <- rep(x, ncur)
		all.idx[[x]] <- seq_len(ncur)
		num.pairs[[x]] <- ncur
	}
	all.a <- unlist(all.a)
	all.t <- unlist(all.t)
	all.mode <- unlist(all.mode)
	all.idx <- unlist(all.idx)
	num.pairs <- unlist(num.pairs)

	# Ordering by anchor, target.
	o <- order(all.a, all.t)
	all.a <- all.a[o]
	all.t <- all.t[o]
	all.mode <- all.mode[o]
	all.idx <- all.idx[o]

	is.diff <- c(TRUE, diff(all.a)!=0L | diff(all.t)!=0L)
	now.index <- cumsum(is.diff)
	by.mode <- split(seq_along(is.diff), all.mode)
	
	indices <- list()
	for (x in seq.it.nk) {
		chosen <- by.mode[[as.character(x)]]
		current.out <- integer(length(chosen))
		current.out[all.idx[chosen]] <- now.index[chosen]
		indices[[x]] <- current.out
	}
	names(indices) <- names(all.hits)
	
	# Selecting the boundaries to report.
	if (minbox) {
		a.chrs <- a.starts <- a.ends <- t.chrs <- t.starts <- t.ends <- list()
		for (x in seq.it.nk) {
			current <- all.hits[[x]]
			aid <- anchors(current, id=TRUE)
			tid <- targets(current, id=TRUE)
			rstarts <- start(regions(current))
			rends <- end(regions(current))
			rchrs <- as.character(seqnames(regions(current)))
			a.chrs[[x]] <- rchrs[aid]
			a.starts[[x]] <- rstarts[aid]
			a.ends[[x]] <- rends[aid] + 1L
			t.chrs[[x]] <- rchrs[tid]
			t.starts[[x]] <- rstarts[tid]
			t.ends[[x]] <- rends[tid] + 1L
		}
		boxed <- .minBoundingBox(unlist(indices), unlist(a.chrs), unlist(a.starts), unlist(a.ends), 
			unlist(t.chrs), unlist(t.starts), unlist(t.ends), seqinfo(parents))
		anchors <- boxed$anchors
		targets <- boxed$targets
	} else {
		anchors <- parents[all.a[is.diff]]
		targets <- parents[all.t[is.diff]]
	}
	return(list(indices=indices, anchors=anchors, targets=targets))
}
