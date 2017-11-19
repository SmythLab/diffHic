###################################################################################################
# This tests the interaction counting capabilities of the patch counter.

Sys.setlocale(category="LC_COLLATE",locale="C")
chromos<-c(chrA=51, chrB=31)
source("simcounts.R")

dir.create("temp-patch")
dir1<-"temp-patch/1.h5"

comp <- function(npairs1, dist, cuts, cap=NA) {
    simgen(dir1, npairs1, chromos)
    param <- pairParam(fragments=cuts, cap=cap)
    y <- squareCounts(dir1, param=param, width=dist, filter=1L)

    # Checking overlaps with one region.
    dummy.1 <- suppressWarnings(resize(regions(y)[sample(length(regions(y)), 1)], fix="center", dist*2))
    ysub <- extractPatch(dir1, param, dummy.1, width=dist)
    m <- overlapsAny(y, dummy.1, use.region="first") & overlapsAny(y, dummy.1, use.region="second")
    ref <- y[m,]

    ref <- sort(ref)
    ysub <- sort(ysub)
    stopifnot(identical(assay(ref), assay(ysub)))
    stopifnot(identical(anchors(ref, id=TRUE), anchors(ysub, id=TRUE)))
    stopifnot(identical(regions(ref), regions(ysub)))

    # Checking overlaps with multiple regions.
    dummy.2 <- suppressWarnings(resize(regions(y)[sample(length(regions(y)), 1)], fix="center", dist*2))
    ysub <- extractPatch(dir1, param, dummy.1, dummy.2, width=dist)
    m <- linkOverlaps(y, dummy.1, dummy.2)
    ref <- y[m$query,]

    ref <- sort(ref)
    ysub <- sort(ysub)
    stopifnot(identical(assay(ref), assay(ysub)))
    stopifnot(identical(anchors(ref, id=TRUE), anchors(ysub, id=TRUE)))
    stopifnot(identical(regions(ref), regions(ysub)))

    # Checking we get the same results with restrict.regions=TRUE.
    y.alt <- extractPatch(dir1, param, dummy.1, dummy.2, width=dist, restrict.regions=TRUE)
    if (!identical(assay(ysub), assay(y.alt)) ||  !identical(anchors(ysub), anchors(y.alt)) ||
            any(! seqlevelsInUse(regions(y.alt)) %in% as.character(c(seqnames(dummy.1), seqnames(dummy.2))))) {
        stop("restrict.regions=TRUE doesn't work for extractPatch")        
    }
    
    output <- interactions(ref)
    output$count <- assay(ref)
    return(output)
}

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

###################################################################################################

set.seed(200)
comp(20, 10000, cuts=simcuts(chromos), cap=NA)
comp(20, 5000, cuts=simcuts(chromos), cap=NA)
comp(20, 20000, cuts=simcuts(chromos), cap=NA)

comp(50, 10000, cuts=simcuts(chromos), cap=NA)
comp(50, 5000, cuts=simcuts(chromos), cap=NA)
comp(50, 20000, cuts=simcuts(chromos), cap=NA)

comp(100, 10000, cuts=simcuts(chromos), cap=NA)
comp(100, 5000, cuts=simcuts(chromos), cap=NA)
comp(100, 20000, cuts=simcuts(chromos), cap=NA)

comp(100, 10000, cuts=simcuts(chromos), cap=1)
comp(100, 5000, cuts=simcuts(chromos), cap=2)
comp(100, 20000, cuts=simcuts(chromos), cap=5)

###################################################################################################
# Wrapping up.

unlink("temp-patch", recursive=TRUE)

###################################################################################################
# End.

