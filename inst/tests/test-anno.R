####################################################################################################
# This tests the clusterPairs function.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

####################################################################################################

simgen <- function(alln, chromos, width, min.space, max.space) {
	# Randomly sampling chromosomes to generate both a pair matrix and the chromosome space.
	output <- GRanges()
	for (x in names(chromos)) {
		n <- chromos[[x]]
		gaps <- round(runif(n, min.space, max.space))
		starts <- cumsum(gaps)
		ends <- starts + width
		suppressWarnings(output <- c(output, GRanges(x, IRanges(starts, ends))))
	}
	
	# Randomly sampling pairs.
   	total <- length(output)
   	chosen1 <- round(runif(alln, 1, total))
   	chosen2 <- round(runif(alln, 1, total))
   	chosen.a <- pmax(chosen1, chosen2)
   	chosen.t <- pmin(chosen1, chosen2)

	# Enforcing uniqueness.
   	o <- order(chosen.a, chosen.t)
   	chosen.a <- chosen.a[o]
   	chosen.t <- chosen.t[o]
   	is.diff <- c(TRUE, diff(chosen.a)!=0 | diff(chosen.t)!=0)
	return(InteractionSet(list(counts=matrix(0L, nrow=alln, ncol=1)), 
        GInteractions(anchor1=chosen.a, anchor2=chosen.t, region=output, mode="reverse"),
        colData=DataFrame(totals=100)))
}

simregs <- function(data) {
    total.size <- range(regions(data))
    out <- tile(range(regions(data)), length(LETTERS))
    for (i in seq_along(out)) {
        names(out[[i]]) <- paste0(LETTERS, i)
    }
    unlist(out)
}

annoref <- function(data, indices, regions, ...) {
    anno1 <- findOverlaps(anchors(data, type="first"), regions, ...)
    anno2 <- findOverlaps(anchors(data, type="second"), regions, ...)
    a1.anno <- sapply(split(names(regions)[subjectHits(anno1)], indices[queryHits(anno1)]), diffHic:::.uniqConcat)
    a2.anno <- sapply(split(names(regions)[subjectHits(anno2)], indices[queryHits(anno2)]), diffHic:::.uniqConcat)
    anchor1.anno <- anchor2.anno <- character(max(indices))
    anchor1.anno[as.integer(names(a1.anno))] <- a1.anno
    anchor2.anno[as.integer(names(a2.anno))] <- a2.anno
    return(list(anchor1=anchor1.anno, anchor2=anchor2.anno))
}

.check_strings <- function(alpha, bravo) {
    if (length(alpha)!=length(bravo)) { stop("strings are not identical") }
    for (x in seq_along(alpha)) { 
        a1 <- sort(unlist(strsplit(alpha[x], split=",")))
        b1 <- sort(unlist(strsplit(bravo[x], split=",")))
        if (!identical(a1, b1)) { stop("strings are not identical") }
    }
    return(NULL)
}

annocomp <- function(data, simregs, indices, ...) {
    # Standard.
    ref <- annoref(data, indices, simregs, ...)
    obs <- annotatePairs(data, indices=indices, regions=simregs, ...)
    .check_strings(ref$anchor1, obs$anchor1)
    .check_strings(ref$anchor2, obs$anchor2)

    # Missing indices.
    ref.x <- annoref(data, seq_along(data), simregs, ...)
    obs.x <- annotatePairs(data, regions=simregs, ...)
    .check_strings(ref.x$anchor1, obs.x$anchor1)
    .check_strings(ref.x$anchor2, obs.x$anchor2)

    # Everyone split up into chunks.
    grouping <- sample(3, length(data), replace=TRUE)
    re.data <- split(data, grouping)
    re.index <- split(indices, grouping)
    obs.i <- annotatePairs(re.data, indices = re.index, regions = simregs, ...)
    .check_strings(ref$anchor1, obs.i$anchor1)
    .check_strings(ref$anchor2, obs.i$anchor2)

    return(head(data.frame(obs)))
}

####################################################################################################

set.seed(3413094)

chromos <- c(chrA=10, chrB=20, chrC=40)
data <- simgen(100, chromos, 20, 50, 100)
regions <- simregs(data)
indices <- sample(5, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

indices <- sample(20, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

indices <- sample(50, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

# Repeating with some different settings.
data <- simgen(100, chromos, 20, 10, 20)
regions <- simregs(data)
indices <- sample(5, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

indices <- sample(20, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

indices <- sample(50, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

# Repeating again, with some different settings.
data <- simgen(100, chromos, 20, 100, 200)
regions <- simregs(data)
indices <- sample(5, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

indices <- sample(20, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

indices <- sample(50, length(data), replace=TRUE)
annocomp(data, regions, indices)
annocomp(data, regions, indices, type="within")

####################################################################################################
# End.
