plotPlaid <- function(file, param, first.region, second.region=first.region,
   	width=10000, col="black", max.count=20, xlab=NULL, ylab=NULL, 
	diag=TRUE, count=FALSE, count.args=list(), ...)
# This function takes a set of boundaries and a count directory and it generates a plaid plot of 
# the result. The plot is colour coded with heat for the desired count range (white for nothing
# and shades of 'col' according to the requested 'alpha'). Binning of restriction fragments is used
# to compute colours, using a bin width of approximately 'width'. Using overlapping rectangles
# should be avoided due to rounding imprecision when computing colours.
#
# written by Aaron Lun
# sometime in 2012.
# last modified 13 May 2017 
{
	first.chr <- as.character(seqnames(first.region))
	second.chr <- as.character(seqnames(second.region))
	if (length(first.chr)!=1L) { stop("exactly one first range is required for plotting") }
	if (length(second.chr)!=1L) { stop("exactly one second range is required for plotting") }

    fragments <- param$fragments
    first.min <- max(1L, start(first.region))
    first.max <- min(seqlengths(fragments)[[first.chr]], end(first.region))
    second.min <- max(1L, start(second.region))
    second.max <- min(seqlengths(fragments)[[second.chr]], end(second.region))
    if (first.min > first.max || second.min > second.max) { stop("invalid anchor ranges supplied") }

    # Stetching a little to allow for some space beyond the plot boundaries.
    f.expanded <- suppressWarnings(trim(resize(first.region, fix="center", width=width(first.region)*2 + 200)))
    s.expanded <- suppressWarnings(trim(resize(second.region, fix="center", width=width(second.region)*2 + 200)))
    patch <- extractPatch(file, param, first.region=first.region, second.region=second.region, width=width)

	# Generating a plot.
	if (is.null(xlab)) { xlab <- first.chr }
	if (is.null(ylab)) { ylab <- second.chr }
	plot(-1, -1, xlim=c(first.min, first.max), ylim=c(second.min, second.max), xlab=xlab, ylab=ylab, type="n", ...)

	# Setting up the color function.		
	my.col<-col2rgb(col)[,1]
	colfun <- function(count) { .get.new.col(my.col, pmin(1, count/max.count)) }
    if (nrow(patch)==0L) { return(invisible(colfun)) }
	
    # Assigning coordinates to x and y-axes, while checking whether chromosome names have been flipped.
	if (!metadata(patch)$flipped) { 
		x.ranges <- anchors(patch, type="first")
		y.ranges <- anchors(patch, type="second")
	} else {
		x.ranges <- anchors(patch, type="second")  
        y.ranges <- anchors(patch, type="first")
	}

	# Summoning a function to get colours.
	all.cols <- colfun(assay(patch))
	labels <- NULL
	if (count) { labels <- assay(patch) } 
	.plotDiag(x.ranges, y.ranges, all.cols, diag=diag, labels=labels, label.args=count.args)
	return(invisible(colfun))
}

###########################################################
# Helper functions.

.plotDiag <- function(xranges, yranges, colors, diag=FALSE, do.label=FALSE, labels=NULL, label.args=list()) {
	for (it in seq_len(2)) {
		# Adding boxes (and text in those boxes, if desired).
		rect(xleft=start(xranges)-0.5, xright=end(xranges)+0.5, ybottom=start(yranges)-0.5, ytop=end(yranges)+0.5, border=NA, col=colors)
		if (!is.null(labels)) { do.call(text, c(list(x=mid(ranges(xranges)), y=mid(ranges(yranges)), labels=labels), label.args)) }

		# If we want to show elements past the diagonal for intra-chromosomal plots, we do so.
		if (!diag || as.character(seqnames(xranges[1]))!=as.character(seqnames(yranges[1]))) { break }
		offdiag <- xranges!=yranges
		if (!any(offdiag)) { break }
		temp <- xranges[offdiag]
		xranges <- yranges[offdiag]
		yranges <- temp
		colors <- colors[offdiag]
	}
	invisible(NULL)
}

.get.new.col <- function(old.colour, strength) {
	if (!length(strength)) { return(character(0)) }
	old.colour <- as.integer(old.colour)
	remnant <- 255 - old.colour
	adj <- outer(remnant, 1-strength) + old.colour
	rgb(adj[1,], adj[2,], adj[3,], maxColorValue=255)
}

###########################################################

plotDI <- function(data, fc, first.region, second.region=first.region,
	col.up="red", col.down="blue",background="grey70", zlim=NULL, 
	xlab=NULL, ylab=NULL, diag=TRUE, ...)
# This function plots differential interactions.
#
# written by Aaron Lun
# created 21 November 2014
# last modified 13 May 2017
{
    # Checking for proper type.
    .check_StrictGI(data)
	
    first.chr <- as.character(seqnames(first.region))
	second.chr <- as.character(seqnames(second.region))
	if (length(first.chr)!=1L) { stop("exactly one first range is required for plotting") }
	if (length(second.chr)!=1L) { stop("exactly one second range is required for plotting") }

	# Setting up the boundaries.
	first.min <- max(1L, start(first.region))
	first.max <- min(seqlengths(regions(data))[[first.chr]], end(first.region))
	second.min <- max(1L, start(second.region))
	second.max <- min(seqlengths(regions(data))[[second.chr]], end(second.region))
	if (first.min > first.max || second.min > second.max) { stop("invalid anchor ranges supplied") }

	# Checking that our points are consistent.
	nr <- nrow(data)
	if (nr!=length(fc)) { stop("length of fold-change vector should equal number of bin pairs") }

	# Identifying the region pairs in our ranges of interest (some stretch, to allow for partial overlaps at plot edge).
	first.keep <- overlapsAny(regions(data), first.region, maxgap=width(first.region)/2+100)
	anchor1.id <- anchors(data, type="first", id=TRUE)
    anchor2.id <- anchors(data, type="second", id=TRUE)
	flipped <- FALSE

	if (suppressWarnings(first.region!=second.region)) {
		second.keep <- overlapsAny(regions(data), second.region, maxgap=width(second.region)/2+100)
		keep.normal <- second.keep[anchor2.id] & first.keep[anchor1.id]
		keep.flipped <- first.keep[anchor2.id] & second.keep[anchor1.id]
		if (!any(keep.normal) && any(keep.flipped)) { flipped <- TRUE } # Checking whether chromosome names are flipped.
		keep <- keep.normal | keep.flipped # Include all bin pairs while accounting for reflection around diagonal.
	} else {
		keep <- first.keep[anchor1.id] & first.keep[anchor2.id]
	}

	if (is.null(xlab)) { xlab <- first.chr }
	if (is.null(ylab)) { ylab <- second.chr }
	plot(-1, -1, xlim=c(first.min, first.max), ylim=c(second.min, second.max), xlab=xlab, ylab=ylab, type="n", ...)
	u <- par("usr") # The coordinates of the plot area
	rect(u[1], u[3], u[2], u[4], col=background, border=NA)
	box()

	# Generating a colour-generating function.
	if (is.null(zlim)) { zlim <- max(abs(fc)) }
	up.col <- col2rgb(col.up)[,1]
	down.col <- col2rgb(col.down)[,1]
	colfun <- .forgeNewFun(zlim, up.col, down.col)
	if (!any(keep))	{ return(invisible(colfun)) }
	
	# Assigning coordinates while checking if the chromosome names are flipped.		
	current <- data[keep,]
	if (flipped) {
		first.ranges <- anchors(current, type="second")
		second.ranges <- anchors(current, type="first")
	} else {
		first.ranges <- anchors(current, type="first")
		second.ranges <- anchors(current, type="second")
	}	

	# Making the actual plot.
	all.cols <- colfun(fc[keep])
	.plotDiag(first.ranges, second.ranges, all.cols, diag=diag)
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
