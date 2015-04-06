getDistance <- function(data, type=c("mid", "gap", "span")) 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
#
# written by Aaron Lun
# created 22 April 2014
# last modified 24 March 2015
{
	aid <- anchors(data, id=TRUE)
	tid <- targets(data, id=TRUE)
	st <- start(regions(data))
	en <- end(regions(data))
	chr <- as.character(seqnames(regions(data)))

	is.same <- chr[aid]==chr[tid]
	all.as <- st[aid[is.same]]
	all.ae <- en[aid[is.same]]
	all.ts <- st[tid[is.same]]
	all.te <- en[tid[is.same]]

	output <- rep(NA, nrow(data))
	type <- match.arg(type)
	if (type=="gap") {
		output[is.same] <- pmax(all.as, all.ts) - pmin(all.ae, all.te) - 1L
	} else if (type=="span") {
		output[is.same] <- pmax(all.ae, all.te) - pmin(all.as, all.ts) + 1L
	} else if (type=="mid") {
		output[is.same] <- as.integer(floor((all.as + all.ae)/2) - floor((all.ts + all.te)/2))
	}
	return(output)
}

getArea <- function(data, bp=TRUE)
# Computes the area of the interaction space, either in base pair-squared terms
# or in terms of pairs of restriction fragments. This allows adjustment of
# abundances for comparison between differently-sized areas. Special behaviour
# is necessary on the diagonal, as reflection halves the space. Coercion to 
# double is necessary to prevent overflows of the integer type.
# 
# written by Aaron Lun
# created 30 July 2014
# last modified 20 March 2015
{
	ax <- anchors(data, id=TRUE)
	tx <- targets(data, id=TRUE)
	reg <- regions(data)

	if (bp) {
		cur.width <- as.double(width(reg))
		returned <- cur.width[ax] * cur.width[tx]

		# Accounting for special behaviour around the diagonal.	It you don't halve,
		# you'll have to double every other (unreflected) area.
		overlap <- getDistance(data, type="gap")
		is.olap <- !is.na(overlap) & overlap < -0.5
		lap.dist <- -overlap[is.olap]
		self.lap.area <- lap.dist * (lap.dist - 1)/2		
		returned[is.olap] <- returned[is.olap] - self.lap.area
	} else {
		is.same <- ax==tx
		curnfrag <- as.double(reg$nfrags[ax])
		returned <- curnfrag * reg$nfrags[tx]
		returned[is.same] <- curnfrag[is.same]*(curnfrag[is.same]+1)/2

		# Detour to protect against overlapping regions.
		left.edge <- pmax(start(reg)[ax], start(reg)[tx])
		right.edge <- pmin(end(reg)[ax], end(reg)[tx])
		is.partial <- !is.same & right.edge >= left.edge & 
			as.logical(seqnames(reg)[ax]==seqnames(reg)[tx]) 

		if (any(is.partial)) { 
			right.edge <- right.edge[is.partial]
			left.edge <- left.edge[is.partial]
			by.chr <- split(1:sum(is.partial), as.character(seqnames(reg)[ax][is.partial]))
			fragments <- exptData(data)$param$fragments
			fdata <- .delimitFragments(fragments)

			for (x in 1:length(fdata$chr)) {
				current.chr <- fdata$chr[x]
				curdex <- by.chr[[current.chr]]
				if (is.null(curdex)) { next }
		
				indices <- fdata$start[x]:fdata$end[x]
				right.olap <- match(right.edge[curdex], end(fragments)[indices])
				left.olap <- match(left.edge[curdex], start(fragments)[indices])
 	    		if (any(is.na(right.olap)) || any(is.na(left.olap))) { stop("region boundaries should correspond to restriction fragment boundaries") }
		
				n.overlap <- right.olap - left.olap + 1	
				returned[is.partial][curdex] <- returned[is.partial][curdex] - n.overlap*(n.overlap-1)/2
			}
		}
	}

	return(returned)
}
