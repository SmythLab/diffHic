rotPlaid <- function(file, param, region, width=10000, col="red", max.count=20, xlab=NULL, ylab="Gap", ...)
# This constructs a sideways plot of interaction intensities.
# Boxes represent interactions where the interacting loci are
# on the x-axis, extended from the diagonal.
#
# written by Aaron Lun
# created 18 September 2014
# last modified 20 March 2015
{
	if (!is.integer(width)) { width<-as.integer(width) }
	xchr <- as.character(seqnames(region))
	xstart <- start(region)
	xend <- end(region)
	if (length(xchr)!=1L) { stop("exactly one region is required for plotting") }

	# Setting up the parameters
	fragments <- param$fragments
	if (!xchr %in% seqlevels(fragments)) { stop("anchor/target chromosome names not in cut site list") } 
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Setting up the boundaries.
	x.min <- max(1L, xstart)
	x.max <- min(seqlengths(fragments)[[xchr]], xend)
	if (x.min >= x.max) { stop("invalid anchor/target ranges supplied") }
						
	# Identifying the fragments in our ranges of interest (with some leeway, to ensure that edges of the plot are retained).
	keep <- overlapsAny(fragments, region, maxgap=width(region)/2)
	new.pts <- .getBinID(fragments[keep], width)
	out.id <- integer(length(fragments))
	out.id[keep] <- new.pts$id
	
	# Pulling out the read pair indices from each file.
	all.dex <- .loadIndices(file, seqlevels(param$fragments))
	if (!is.null(all.dex[[xchr]][[xchr]])) {
		current <- .baseHiCParser(TRUE, file, xchr, xchr, discard=discard, cap=cap)[[1]]
	} else { 
		current<-data.frame(anchor.id=integer(0), target.id=integer(0))
	}

	# Computing the max height.
	max.height <- x.max - x.min
	if (is.null(xlab)) { xlab <- xchr }
	plot(-1, -1, xlim=c(x.min, x.max), ylim=c(0, max.height),
		xlab=xlab, yaxs="i", ylab=ylab, type="n", bg="transparent", ...)
	if (!nrow(current))	{ next }

   	retain <- keep[current$anchor.id] & keep[current$target.id]
	out<-.Call(cxx_count_patch, list(current[retain,]), out.id, 1L)
	if (is.character(out)) { stop(out) }

	# Rotating the vertices.
	anchors <- new.pts$region[out[[1]]]
	targets <- new.pts$region[out[[2]]]
	corner <- .spawnVertices(anchors, targets)

	# Plotting these new vertices.
	my.col<-col2rgb(col)[,1]
	colfun <- function(count) { .get.new.col(my.col, pmin(1, count/max.count)) }
	polygon(corner$x, corner$y, border=NA, col=colfun(out[[3]]))
	return(invisible(colfun))
}

#################################################################

.spawnVertices <- function(anchors, targets) {
	all.x <- all.y <- rep(NA, length(anchors)*5L-1L)
	hits <- 0:(length(anchors)-1L) * 5L
	counter <- 1L
	for (mode in list(c(1,1), c(1,2), c(2,2), c(2,1))) {
		if (mode[1]==1L) { 
			cur.x <- start(anchors) - 0.5
		} else {
			cur.x <- end(anchors) + 0.5
		} 
		if (mode[2]==1L) { 
			cur.y <- start(targets) - 0.5
		} else {
			cur.y <- end(targets) + 0.5
		}
		all.x[counter+hits] <- (cur.x + cur.y)/2
		all.y[counter+hits] <- cur.x - cur.y
		counter <- counter + 1L
	}
	return(list(x=all.x, y=all.y))
}

#################################################################

rotDI <- function(data, fc, region, col.up="red", col.down="blue",
    background="grey70", zlim=NULL, xlab=NULL, ylab="Gap", ...)
# This constructs a sideways plot of interaction intensities.
# Boxes represent interactions where the interacting loci are
# on the x-axis, extended from the diagonal.
#
# written by Aaron Lun
# created 18 September 2014
# last modified 20 March 2015
{
	xchr <- as.character(seqnames(region))
	xstart <- start(region)
	xend <- end(region)
	if (length(xchr)!=1L) { stop("exactly one region is required for plotting") }

	# Setting up the boundaries.
	x.min <- max(1L, xstart)
	x.max <- min(seqlengths(regions(data))[[xchr]], xend)
	if (x.min >= x.max) { stop("invalid anchor/target ranges supplied") }
						
	# Identifying the fragments in our ranges of interest (with some leeway, to ensure that edges of the plot are retained).
	ref.keep <- overlapsAny(regions(data), region, maxgap=width(region)/2)
	keep <- ref.keep[anchors(data, id=TRUE)] & ref.keep[targets(data, id=TRUE)]

	# Computing the max height.
	max.height <- x.max - x.min
	if (is.null(xlab)) { xlab <- xchr }
	plot(-1, -1, xlim=c(x.min, x.max), ylim=c(0, max.height), xlab=xlab, yaxs="i", ylab=ylab, type="n", ...)
	u <- par("usr") # The coordinates of the plot area
	rect(u[1], u[3], u[2], u[4], col=background, border=NA)
	if (!any(keep)) { next }

	# Generating colours.
	if (is.null(zlim)) { zlim <- max(abs(fc)) }
	up.col <- col2rgb(col.up)[,1]
	down.col <- col2rgb(col.down)[,1]
	colfun <- .forgeNewFun(zlim, up.col, down.col)

	# Rotating the vertices and plotting them..
	current <- data[keep,]
	corner <- .spawnVertices(anchors(current), targets(current))
	polygon(corner$x, corner$y, border=NA, col=colfun(fc[keep]))
	return(invisible(colfun))
}

