rotPlaid <- function(file, param, region, width=10000, col="black", max.count=20, xlab=NULL, max.height=NULL, ylab="Gap", ...)
# This constructs a sideways plot of interaction intensities.
# Boxes represent interactions where the interacting loci are
# on the x-axis, extended from the diagonal.
#
# written by Aaron Lun
# created 18 September 2014
# last modified 22 November 2015
{
	xchr <- as.character(seqnames(region))
	xstart <- start(region)
	xend <- end(region)
	if (length(xchr)!=1L) { stop("exactly one region is required for plotting") }

	# Setting up the parameters
	fragments <- param$fragments
	if (!xchr %in% seqlevelsInUse(fragments)) { stop("chromosome name not in cut site list") } 
	discard <- .splitDiscards(param$discard)
	cap <- param$cap
	frag.by.chr <- .splitByChr(fragments)

	# Setting up the boundaries.
	x.min <- max(1L, xstart)
	x.max <- min(seqlengths(fragments)[[xchr]], xend)
	if (x.min >= x.max) { stop("invalid ranges supplied") }
    if (is.null(max.height)) { max.height <- x.max - x.min }

	# Setting up the boxes.		
	width<-as.integer(width) 
	cur.chrs <- frag.by.chr$first[[xchr]]:frag.by.chr$last[[xchr]]
	new.pts <- .getBinID(fragments[cur.chrs], width)
	out.id <- integer(length(fragments))
	out.id[cur.chrs] <- new.pts$id
						
	# Identifying the boxes in our ranges of interest (with some leeway, to ensure that 
	# there's stuff in the corners of the rotated plot). Specifically, you need to include 
	# 'center +/- max.height' on either side to fill up the top left/right corners; this is
	# equivalent to the region interval +- 'max.height/2'. We ask for a bit more, to be safe.
	use.bin <- overlapsAny(new.pts$region, region, maxgap=max.height*0.7)
	keep.frag <- logical(length(fragments))
	keep.frag[cur.chrs] <- use.bin[new.pts$id]	
	
	# Pulling out the read pair indices from each file.
	all.dex <- .loadIndices(file, seqlevelsInUse(fragments))
	if (!is.null(all.dex[[xchr]][[xchr]])) {
		current <- .baseHiCParser(TRUE, file, xchr, xchr, chr.limits=frag.by.chr,
			discard=discard, cap=cap)[[1]]
	} else { 
		current<-data.frame(anchor1.id=integer(0), anchor2.id=integer(0))
	}

	# Making the plot.
	if (is.null(xlab)) { xlab <- xchr }
	plot(-1, -1, xlim=c(x.min, x.max), ylim=c(0, max.height),
		xlab=xlab, yaxs="i", ylab=ylab, type="n", bg="transparent", ...)

	# Getting the colour (and returning it, if necessary).
	my.col<-col2rgb(col)[,1]
	colfun <- function(count) { .get.new.col(my.col, pmin(1, count/max.count)) }

	# Collating read pairs into counts.
   	retain <- keep.frag[current$anchor1.id] & keep.frag[current$anchor2.id]
	if (!any(retain)) { return(invisible(colfun)) }
	bin.indices <- out.id[keep.frag]
	out<-.Call(cxx_count_patch, list(current[retain,]), out.id, 1L, 
			bin.indices[1L], tail(bin.indices, 1L)) # First and last bin indices on interval of interest.
	if (is.character(out)) { stop(out) }

	# Rotating the vertices.
	anchor1.regions <- new.pts$region[out[[1]]]
    anchor2.regions <- new.pts$region[out[[2]]]
	corner <- .spawnVertices(anchor1.regions, anchor2.regions)

	# Plotting these new vertices.
	polygon(corner$x, corner$y, border=NA, col=colfun(out[[3]]))
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
# last modified 8 December 2015
{
    # Checking for proper type.
    .check_StrictGI(data)

	xchr <- as.character(seqnames(region))
	xstart <- start(region)
	xend <- end(region)
	if (length(xchr)!=1L) { stop("exactly one region is required for plotting") }

	# Setting up the boundaries.
	x.min <- max(1L, xstart)
	x.max <- min(seqlengths(regions(data))[[xchr]], xend)
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

