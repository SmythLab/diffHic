rotPlaid <- function(file, param, region, width=10000, col="black", max.count=20, xlab=NULL, max.height=NULL, ylab="Gap", ...)
# This constructs a sideways plot of interaction intensities.
# Boxes represent interactions where the interacting loci are
# on the x-axis, extended from the diagonal.
#
# written by Aaron Lun
# created 18 September 2014
# last modified 17 March 2017
{
	xchr <- as.character(seqnames(region))
	if (length(xchr)!=1L) { stop("exactly one region is required for plotting") }

    fragments <- param$fragments
	x.min <- max(1L, start(region))
	x.max <- min(seqlengths(fragments)[[xchr]], end(region))
	if (x.min >= x.max) { stop("invalid ranges supplied") }
    if (is.null(max.height)) { max.height <- x.max - x.min }

    # Expanding the region to account for top-right/left regions of the plot.
    expanded <- suppressWarnings(trim(resize(region, fix="center", width=width(region) + max.height*1.5)))
    patch <- extractPatch(file, param, first.region=expanded, width=width)

	# Making the plot.
	if (is.null(xlab)) { xlab <- xchr }
	plot(-1, -1, xlim=c(x.min, x.max), ylim=c(0, max.height),
		xlab=xlab, yaxs="i", ylab=ylab, type="n", bg="transparent", ...)

	# Getting the colour (and returning it, if necessary).
	my.col<-col2rgb(col)[,1]
	colfun <- function(count) { .get.new.col(my.col, pmin(1, count/max.count)) }

	# Rotating the vertices.
	anchor1.regions <- anchors(patch, type="first")
    anchor2.regions <- anchors(patch, type="second")
	corner <- .spawnVertices(anchor1.regions, anchor2.regions)

	# Plotting these new vertices.
	polygon(corner$x, corner$y, border=NA, col=colfun(assay(patch)))
	return(invisible(colfun))
}

#################################################################

.spawnVertices <- function(anchor1.regions, anchor2.regions) {
	all.x <- all.y <- rep(NA, length(anchor1.regions)*5L-1L)
	hits <- 0:(length(anchor1.regions)-1L) * 5L
	counter <- 1L
	for (mode in list(c(1,1), c(1,2), c(2,2), c(2,1))) {
		if (mode[1]==1L) { 
			cur.x <- start(anchor1.regions) - 0.5
		} else {
			cur.x <- end(anchor1.regions) + 0.5
		} 
		if (mode[2]==1L) { 
			cur.y <- start(anchor2.regions) - 0.5
		} else {
			cur.y <- end(anchor2.regions) + 0.5
		}
		all.x[counter+hits] <- (cur.x + cur.y)/2
		all.y[counter+hits] <- cur.x - cur.y
		counter <- counter + 1L
	}
	return(list(x=all.x, y=all.y))
}

#################################################################

rotDI <- function(data, fc, region, col.up="red", col.down="blue",
	background="grey70", zlim=NULL, xlab=NULL, max.height=NULL, ylab="Gap", ...)
# This constructs a sideways plot of interaction intensities.
# Boxes represent interactions where the interacting loci are
# on the x-axis, extended from the diagonal.
#
# written by Aaron Lun
# created 18 September 2014
# last modified 13 May 2017
{
    # Checking for proper type.
    .check_StrictGI(data)

	# Setting up the boundaries.
	xchr <- as.character(seqnames(region))
	if (length(xchr)!=1L) { stop("exactly one region is required for plotting") }
	x.min <- max(1L, start(region))
	x.max <- min(seqlengths(regions(data))[[xchr]], end(region))
	if (x.min >= x.max) { stop("invalid ranges supplied") }
    if (is.null(max.height)) { max.height <- x.max - x.min }

    # Checking that our points are consistent.
    nr <- nrow(data)
    if (nr!=length(fc)) { stop("length of fold-change vector should equal number of bin pairs") }

    # Identifying the fragments in our ranges of interest.
	ref.keep <- overlapsAny(regions(data), region, maxgap=max.height*0.7)
	keep <- ref.keep[anchors(data, type="first", id=TRUE)] & ref.keep[anchors(data, type="second", id=TRUE)]

	# Making the plot.
	if (is.null(xlab)) { xlab <- xchr }
	plot(-1, -1, xlim=c(x.min, x.max), ylim=c(0, max.height), xlab=xlab, yaxs="i", ylab=ylab, type="n", ...)
	u <- par("usr") # The coordinates of the plot area
	rect(u[1], u[3], u[2], u[4], col=background, border=NA)
	box()

	# Generating colours.
	if (is.null(zlim)) { zlim <- max(abs(fc)) }
	up.col <- col2rgb(col.up)[,1]
	down.col <- col2rgb(col.down)[,1]
	colfun <- .forgeNewFun(zlim, up.col, down.col)
	if (!any(keep)) { return(invisible(colfun)) }

	# Rotating the vertices and plotting them.
	current <- data[keep,]
	corner <- .spawnVertices(anchors(current, type="first"), anchors(current, type="second"))
	polygon(corner$x, corner$y, border=NA, col=colfun(fc[keep]))
	return(invisible(colfun))
}

