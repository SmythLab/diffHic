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

