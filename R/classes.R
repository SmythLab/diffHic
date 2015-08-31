# Defines the DIList class, that will be used to hold various information.

setClass("DIList", representation(counts="matrix", 
			anchors="integer", targets="integer", regions="GRanges",
			colData="DataFrame", exptData="List"))

setValidity("DIList", function(object) {
	if (nrow(object@counts)!=length(object@anchors)) {
		return('rows in count matrix not equal to length of anchor vector')
	} 
	if (nrow(object@counts)!=length(object@targets)) { 
		return('rows in count matrix not equal to length of target vector')
	}
	if (ncol(object@counts)!=nrow(object@colData)) {
		return('columns of count matrix not equal to rows of column data frame')
	}

	if (!all(object@anchors >= 1L)) { 
		return('not all anchors are positive integers')
	} 
	if (!all(object@targets >= 1L)) {
		return('not all targets are positive integers')
	}
	if (!all(object@anchors <= length(object@regions))) {
		return('not all anchors refer to valid regions')
	} 
	if (!all(object@targets <= length(object@regions))) { 
		return('not all targets refer to valid regions')
	}
	if (!all(object@anchors >= object@targets)) { 
		return('target indices cannot be greater than anchor indices')
	}
	return(TRUE)
})

setMethod("initialize", signature("DIList"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("show", signature("DIList"), function(object) {
	total <- nrow(object@counts)
	nregs <- length(object@regions)
	nlibs <- ncol(object@counts)
	cat("DIList object for", nlibs, ifelse(nlibs==1L, "library", "libraries"), 
		"with", total, ifelse(total==1L, "pair", "pairs"), "across", 
		nregs, ifelse(nregs==1L, "region\n", "regions\n"))
})

# Assorted subsetting methods.
setMethod("[", "DIList", function(x, i, j, ..., drop=TRUE) {
	if (missing(i)) {
		new.counts <- x@counts				
		new.anchors <- x@anchors
		new.targets <- x@targets	
	} else {
		new.counts <- x@counts[i,,drop=FALSE]	
		new.anchors <- x@anchors[i]
		new.targets <- x@targets[i]
	}

	if (missing(j)) { 
		new.colData <- x@colData
	} else {
		new.counts <- new.counts[,j,drop=FALSE]
		new.colData <- x@colData[j,,drop=FALSE]
	}
	initialize(x, counts=new.counts, colData=new.colData, 
		anchors=new.anchors, targets=new.targets,
		regions=x@regions, exptData=x@exptData)
})

# Some getters. No need for setters, really.
setGeneric("anchors", function(object, ...) { standardGeneric("anchors") })
setMethod("anchors", signature("DIList"), function(object, id=FALSE) {
	if (id) { return(object@anchors) }
	object@regions[object@anchors]
})

setGeneric("targets", function(object, ...) { standardGeneric("targets") })
setMethod("targets", signature("DIList"), function(object, id=FALSE) {
	if (id) { return(object@targets) }
	object@regions[object@targets]
})

#setGeneric("counts", function(object) { standardGeneric("counts") })
setMethod("counts", signature("DIList"), function(object) {
	object@counts
})

setGeneric("regions", function(object) { standardGeneric("regions") })
setMethod("regions", signature("DIList"), function(object) {
	object@regions
})

setMethod("dim", signature("DIList"), function(x) {
	dim(x@counts)
})

setMethod("dimnames", signature("DIList"), function(x) {
	dimnames(x@counts)
})

setMethod("$", signature("DIList"), function(x, name) { 
	x@colData[[name]]
})

# Borrowing these from GenomicRanges.
setMethod("colData", signature("DIList"), function(x, ...) {
	x@colData
})

setMethod("exptData", signature("DIList"), function(x, ...) {
	x@exptData
})

# Modifier functions.

setMethod("$<-", signature("DIList"), function(x, name, value) { 
	x@colData[[name]] <- value
	x
})

setMethod("exptData<-", signature("DIList", "SimpleList"), function(x, ..., value) {
	x@exptData <- value
	x
})

# Constructor object.
DIList <- function(counts, totals=colSums(counts), anchors, targets, regions, exptData=List(), ...) {
	if (!is.integer(counts)) { storage.mode(counts) <- "integer" }
	anchors <- as.integer(anchors)
	targets <- as.integer(targets)
	totals <- as.integer(totals)
	new("DIList", counts=counts, anchors=anchors, targets=targets, regions=regions,
		colData=DataFrame(totals=totals, ...), exptData=exptData)
}

setMethod("c", signature("DIList"), function (x, ..., add.totals=TRUE, recursive=FALSE) {
	if (!identical(recursive, FALSE)) { 
		stop("recursive argument not supported")
	}
	output <- list(counts(x))
	out.a <- list(anchors(x, id=TRUE))
	out.t <- list(targets(x, id=TRUE))
	totality <- x$totals

	ix <- 2L
	for (i in list(...)) {
		if (!is(i, "DIList")) { 
			stop("elements to be concatenated must be DIList objects")
		}
		if (!identical(regions(x), regions(i))) {
			stop("regions should be identical between DIList objects")
		}

		out.a[[ix]] <- anchors(i, id=TRUE)
		out.t[[ix]] <- targets(i, id=TRUE)
		output[[ix]] <- counts(i)

		if (add.totals) { 
			totality <- totality + i$totals 
		} else if (!identical(i$totals, totality)) { 
			warning("totals are not identical between DIList objects")
		}
		ix <- ix + 1L
	}
	
	colData <- colData(x)
	colData$totals <- totality
	new("DIList", counts=do.call(rbind, output), 
		anchors=unlist(out.a), targets=unlist(out.t), regions=regions(x),
		exptData=exptData(x), colData=colData)
})

setMethod("as.matrix", signature("DIList"), function(x, first=NULL, second=first, fill=NULL, ...) {
	if (!is.null(first)) { keep.first <- as.logical(seqnames(regions(x)) %in% first) }
	else { keep.first <- !logical(length(regions(x))) }
	if (!is.null(second)) { keep.second <- as.logical(seqnames(regions(x)) %in% second) }
	else { keep.second <- !logical(length(regions(x))) }
	new.f <- cumsum(keep.first)
	new.s <- cumsum(keep.second)

	mat <- matrix(NA, nrow=sum(keep.first), ncol=sum(keep.second))
	rownames(mat) <- which(keep.first)
	colnames(mat) <- which(keep.second)
	aid <- anchors(x, id=TRUE)
	tid <- targets(x, id=TRUE)

	retain <- keep.first[aid] & keep.second[tid]
	ax <- new.f[aid[retain]]
	tx <- new.s[tid[retain]]
	flip.retain <- keep.first[tid] & keep.second[aid]
	flip.ax <- new.f[tid[flip.retain]]
	flip.tx <- new.s[aid[flip.retain]]

	if (is.null(fill)) { 
		fill <- numeric(nrow(x))
		retained <- retain | flip.retain
		fill[retained] <- aveLogCPM(asDGEList(x[retained,])) 
	}
	mat[ax + (tx-1L) * nrow(mat)] <- fill[retain] 
	mat[flip.ax + (flip.tx-1L) * nrow(mat)] <- fill[flip.retain] 

	return(mat)
})

# Setting some methods inspired by equivalents in csaw.
setMethod("asDGEList", signature("DIList"), function(object, lib.sizes, ...) {
	if (missing(lib.sizes)) { 
	    if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals 
	}
	DGEList(counts=counts(object), lib.size=lib.sizes, ...)
})

setMethod("normOffsets", signature("DIList"), function(object, lib.sizes, ...) {
	if (missing(lib.sizes)) { 
	    if (is.null(object$totals)) { warning("library sizes not found in 'totals', setting to NULL") }
		lib.sizes <- object$totals 
	}
	normOffsets(counts(object), lib.sizes=lib.sizes, ...)
})

setMethod("normalize", signature("DIList"), function(object, lib.sizes, ...) {
    .Deprecated(new="normOffsets", old="normalize")
	normOffsets(object, lib.sizes=lib.sizes, ...)
})

########################################################################################
# Defining the pairParam class.

setClass("pairParam", representation(fragments="GRanges", restrict="character", discard="GRanges", cap="integer"))

setValidity("pairParam", function(object) {
	# Checking that the fragments are in some order by chromosome name, 
	# and in the correct order by restriction fragment width.		
	if (length(object@fragments)>1L) { 
		if (anyDuplicated(runValue(seqnames(object@fragments)))) { 
			return('restriction fragments should be sorted by chromosome name')	
		}
	
		# <, not <=, to allow nested fragments at start and end of the chromosome when 
		# overhang and site lengths are equal.
		unsort <- diff(start(object@fragments)) < 0L | diff(end(object@fragments)) < 0L 

		# Should be +1, to get to the first element of each chromosome; but, unsort 
		# is missing the first element (because of diff), so no need to add 1.
		unsort[head(cumsum(runLength(seqnames(object@fragments))),-1L)] <- FALSE 

		if (any(unsort)) {
			return('restriction fragments should be sorted by start and end coordinates')
		}
	}

	if (any(strand(object@fragments)!="*") ) {
		return('restriction fragment ranges should be unstranded')
	}

	if (length(object@cap)!=1L || (!is.na(object@cap) && object@cap <= 0L)) { 
		return('any specified cap should be a positive integer')
	}

	if (attributes(object@restrict)$paired && length(object@restrict)%%2!=0L) {
		return('restrict vector must be of even length for paired extraction')
	}
	return(TRUE)
})

setMethod("initialize", signature("pairParam"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("$", signature("pairParam"), function(x, name) { 
	slot(x, name)
})

setMethod("show", signature("pairParam"), function(object) {
#	if (is.na(object@min.inward)) { 
#		cat("No minimum insert size specified for inward-facing read pairs\n")
#	} else {
#		cat("Minimum insert size for inward-facing read pairs is", object@min.inward, "bp\n")
#	} 
#	if (is.na(object@min.outward)) { 
#		cat("No minimum insert size specified for outward-facing read pairs\n")
#	} else {
#		cat("Minimum insert size for outward-facing read pairs is", object@min.outward, "bp\n")
#	}
#	if (is.na(object@max.frag)) {
#		cat("No maximum fragment size specified\n")
#	} else {
#		cat("Maximum fragment size is", object@max.frag, "bp\n")
#	}

	nfrags <- length(object@fragments)
	nchrs <- length(runValue(seqnames(object@fragments)))
	cat("Genome contains", nfrags, "restriction", ifelse(nfrags==1L, "fragment", "fragments"), 
		"across", nchrs, ifelse(nchrs==1L, "chromosome\n", "chromosomes\n"))

	ndisc <- length(object@discard)
	if (!ndisc) { 
		cat("No discard regions are specified\n")
	} else {
		cat(ndisc, ifelse(ndisc==1L, "region", "regions"), "specified in which alignments are discarded\n")
	}

	nr <- length(object@restrict)
	if (!nr) { 
		cat("No limits on chromosomes for read extraction\n")
	} else {
		if (!attributes(object@restrict)$paired) {
			cat("Read extraction is limited to", nr, ifelse(nr==1L, "chromosome\n", "chromosomes\n"))
		} else {
			nr <- length(object@restrict)/2
			seq.it.nr <- seq_len(nr)
			cat("Read extraction is limited to pairs between:\n", 
				paste0("\t'", object@restrict[seq.it.nr], "' and '",
				object@restrict[nr+seq.it.nr], "'\n"), sep="")
		}
	}

	if (is.na(object@cap)) {
 	    cat("No cap on the read pairs per pair of restriction fragments\n")
	} else {
		cat("Cap of", object@cap, "on the read pairs per pair of restriction fragments\n")
	}
})

pairParam <- function(fragments, 
#	min.inward=NA, min.outward=NA, max.frag=NA, 
	discard=GRanges(), restrict=NULL, cap=NA)
# This creates a SimpleList of parameter objects, specifying
# how reads should be extracted from the BAM files. The aim is
# to synchronize read loading throughout the package, such that
# you don't have to manually respecify them in each function.
#
# written by Aaron Lun
# 1 September 2014
{
#	max.frag <- as.integer(max.frag)
#	min.inward <- as.integer(min.inward)
#	min.outward <- as.integer(min.outward)
	restrict <- .editRestrict(restrict) 
	cap <- as.integer(cap)
	new("pairParam", 
#			max.frag=max.frag, min.inward=min.inward, min.outward=min.outward,
		restrict=restrict, discard=discard, fragments=fragments, cap=cap)
}

.editRestrict <- function(restrict) {
	paired <- FALSE
	if (!is.null(dim(restrict))) { 
		if (ncol(restrict)!=2L) { stop("restrict matrix must have two columns") }
		paired <- TRUE
	}
	restrict <- as.character(restrict)
	attr(restrict, "paired") <- paired
	restrict
}

setGeneric("reform", function(x, ...) { standardGeneric("reform") })
setMethod("reform", signature("pairParam"), function(x, ...) {
	incoming <- list(...)
	sn <- slotNames(x)
	for (sx in names(incoming)) {
		val <- incoming[[sx]]
		sx <- match.arg(sx, sn)
		incoming[[sx]] <- switch(sx, 
#			max.frag=as.integer(val),
#			min.inward=as.integer(val),
#			min.outward=as.integer(val),
			restrict=.editRestrict(val),
			cap=as.integer(val),
			val)
	}
	do.call(initialize, c(x, incoming))
}) 

