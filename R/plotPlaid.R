plotPlaid <- function(file, param, anchor, target=anchor, 
   	width=10000, col="red", max.count=20, xlab=NULL, ylab=NULL, 
	diag=TRUE, count=FALSE, count.args=list(), ...)
# This function takes a set of boundaries and a count directory and it generates a plaid plot of 
# the result. The plot is colour coded with heat for the desired count range (white for nothing
# and shades of 'col' according to the requested 'alpha'). Binning of restriction fragments is used
# to compute colours, using a bin width of approximately 'width'. Using overlapping rectangles
# should be avoided due to rounding imprecision when computing colours.
#
# written by Aaron Lun
# sometime in 2012.
# last modified 20 March 2015
{
	width<-as.integer(width) 
	achr <- as.character(seqnames(anchor))
	tchr <- as.character(seqnames(target))
	astart <- start(anchor)
	aend <- end(anchor)
	tstart <- start(target)
	tend <- end(target)
    if (length(achr)!=1L) { stop("exactly one anchor range is required for plotting") }
    if (length(tchr)!=1L) { stop("exactly one target range is required for plotting") }

	# Setting up the parameters
	fragments <- param$fragments
	if (!(achr %in% seqlevels(fragments)) || !(tchr %in% seqlevels(fragments))) { stop("anchor/target chromosome names not in cut site list") }
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Setting up the boundaries.
	a.min <- max(1L, astart)
	a.max <- min(seqlengths(fragments)[[achr]], aend)
	t.min <- max(1L, tstart)
	t.max <- min(seqlengths(fragments)[[tchr]], tend)
	if (a.min >= a.max || t.min >= t.max) { stop("invalid anchor/target ranges supplied") }
	
	# Identifying the fragments in our ranges of interest (with some leeway, to ensure that edges of the plot are retained).
	keep.a <- overlapsAny(fragments, anchor, maxgap=width(anchor)/2)
	if (achr!=tchr || astart!=tstart || aend!=tend) {
		keep.t <- overlapsAny(fragments, target, maxgap=width(target)/2)
		keep <- keep.a | keep.t
	} else {
		keep <- keep.t <- keep.a
	}
	new.pts <- .getBinID(fragments[keep], width)
	out.id <- integer(length(fragments))
	out.id[keep] <- new.pts$id

	# Pulling out the read pair indices from each file.
	all.dex <- .loadIndices(file, seqlevels(param$fragments))
	flipped <- FALSE
	if (!is.null(all.dex[[achr]][[tchr]])) {
		current <- .baseHiCParser(TRUE, file, achr, tchr, discard=discard, cap=cap)[[1]]
	} else if (!is.null(all.dex[[tchr]][[achr]])) { 
		current <- .baseHiCParser(TRUE, file, tchr, achr, discard=discard, cap=cap)[[1]]
		flipped <- TRUE
	} else { current<-data.frame(anchor.id=integer(0), target.id=integer(0)) }

	# Generating a plot.
	if (is.null(xlab)) { xlab <- achr }
	if (is.null(ylab)) { ylab <- tchr }
    plot(-1, -1, xlim=c(a.min, a.max), ylim=c(t.min, t.max), xlab=xlab, ylab=ylab, type="n", ...)
	if (!nrow(current))	{ next }

	# Getting the counts around the area of interest, and then collating them.
   	if (flipped) {
		filter.a <- keep.t
		filter.t <- keep.a
   	} else {
 	   	filter.a <- keep.a
		filter.t <- keep.t	
   	}	   
   	retain <- filter.a[current$anchor.id] & filter.t[current$target.id]
    out<-.Call(cxx_count_patch, list(current[retain,]), out.id, 1L)
    if (is.character(out)) { stop(out) }
	
	# Checking whether it's been flipped, to determine the appropriate coordinates..
	if (flipped) { 
		targets <- new.pts$region[out[[1]]]
		anchors <- new.pts$region[out[[2]]]
	} else {
		anchors <- new.pts$region[out[[1]]]
		targets <- new.pts$region[out[[2]]]
	}

	# Summoning a function to get colours.
	my.col<-col2rgb(col)[,1]
	colfun <- function(count) { .get.new.col(my.col, pmin(1, count/max.count)) }
	all.cols <- colfun(out[[3]])
	labels <- NULL
	if (count) { labels <- out[[3]] }
	.plotDiag(anchors, targets, all.cols, diag=diag, labels=labels, label.args=count.args)
	return(invisible(colfun))
}

###########################################################
# Helper functions.

.plotDiag <- function(anchors, targets, colors, diag=FALSE, do.label=FALSE, labels=NULL, label.args=list()) {
	for (it in 1:2) {
		# Adding boxes (and text in those boxes, if desired).
		rect(xleft=start(anchors)-0.5, xright=end(anchors)+0.5, ybottom=start(targets)-0.5, ytop=end(targets)+0.5, border=NA, col=colors)
		if (!is.null(labels)) { do.call(text, c(list(x=mid(ranges(anchors)), y=mid(ranges(targets)), labels=labels), label.args)) }

		# If we want to show elements past the diagonal for intra-chromosomal plots, we do so.
		if (!diag || as.character(seqnames(anchors[1]))!=as.character(seqnames(targets[1]))) { break }
		offdiag <- anchors!=targets
		if (!any(offdiag)) { break }
		temp <- anchors[offdiag]
		anchors <- targets[offdiag]
		targets <- temp
		colors <- colors[offdiag]
	}
	invisible(NULL)
}

.get.new.col <- function(old.colour, strength) {
	old.colour <- as.integer(old.colour)
	remnant <- 255 - old.colour
	adj <- outer(remnant, 1-strength) + old.colour
	rgb(adj[1,], adj[2,], adj[3,], maxColorValue=255)
}

###########################################################

plotDI <- function(data, fc, anchor, target=anchor, col.up="red", col.down="blue",
 	background="grey70", zlim=NULL, xlab=NULL, ylab=NULL, diag=TRUE, ...)
# This function plots differential interactions.
#
# written by Aaron Lun
# created 21 November 2014
# last modified 2 March 2015
{
	achr <- as.character(seqnames(anchor))
	tchr <- as.character(seqnames(target))
	astart <- start(anchor)
	aend <- end(anchor)
	tstart <- start(target)
	tend <- end(target)
    if (length(achr)!=1L) { stop("exactly one anchor range is required for plotting") }
    if (length(tchr)!=1L) { stop("exactly one target range is required for plotting") }

	# Setting up the boundaries.
	a.min <- max(1L, astart)
	a.max <- min(seqlengths(regions(data))[[achr]], aend)
	t.min <- max(1L, tstart)
	t.max <- min(seqlengths(regions(data))[[tchr]], tend)
	if (a.min >= a.max || t.min >= t.max) { stop("invalid anchor/target ranges supplied") }

	# Checking that our points are consistent.
	nr <- nrow(data)
	if (nr!=length(fc)) { stop("length of fold-change vector should equal number of bin pairs") }

	# Identifying the region pairs in our ranges of interest (some stretch, to allow for partial overlaps).
	a.keep <- overlapsAny(regions(data), anchor, maxgap=width(anchor)/2)
	aid <- anchors(data, id=TRUE)
	tid <- targets(data, id=TRUE)
	if (achr!=tchr || astart!=tstart || aend!=tend) {
		t.keep <- overlapsAny(regions(data), target, maxgap=width(target)/2)
		keep <- (a.keep[aid] & t.keep[tid]) | (a.keep[tid] & t.keep[aid])
	} else {
		keep <- a.keep[aid] & a.keep[tid]
	}

	if (is.null(xlab)) { xlab <- achr }
	if (is.null(ylab)) { ylab <- tchr }
    plot(-1, -1, xlim=c(a.min, a.max), ylim=c(t.min, t.max), xlab=xlab, ylab=ylab, type="n", ...)
	u <- par("usr") # The coordinates of the plot area
	rect(u[1], u[3], u[2], u[4], col=background, border=NA)
	if (!any(keep))	{ next }

	# Generating a colour-generating function.
	if (is.null(zlim)) { zlim <- max(abs(fc)) }
	up.col <- col2rgb(col.up)[,1]
	down.col <- col2rgb(col.down)[,1]
	colfun <- .forgeNewFun(zlim, up.col, down.col)
	
	# Making the actual plot.
	current <- data[keep,]
	all.cols <- colfun(fc[keep])
	.plotDiag(anchors(current), targets(current), all.cols, diag=diag)
	return(invisible(colfun))
}

.forgeNewFun <- function(zlim, up.col, down.col) {
	function(FC) {
		all.cols <- character(length(FC))
		pos.FC <- FC > 0				
		if (any(pos.FC)) { all.cols[pos.FC] <- .get.new.col(up.col, pmin(1, FC[pos.FC]/zlim)) }
		if (!all(pos.FC)) { all.cols[!pos.FC] <- .get.new.col(down.col, pmin(1, -FC[!pos.FC]/zlim)) }
		return(all.cols)
	}
}
