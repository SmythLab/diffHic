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
setMethod("anchors", signature("DIList"), function(x, id=FALSE) {
	if (id) { return(x@anchors) }
	x@regions[x@anchors]
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

setMethod("regions", signature("DIList"), function(x) {
	x@regions
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
	.Deprecated(new="InteractionSet", old="DIList")
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

DI2IS <- function(x) {
    InteractionSet(list(counts=counts(x)), 
                   GInteractions(anchor1=anchors(x, id=TRUE), anchor2=targets(x, id=TRUE), regions=regions(x), mode="reverse"),
                   colData=colData(x), metadata=exptData(x))
}
