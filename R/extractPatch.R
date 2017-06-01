extractPatch <- function(file, param, first.region, second.region=first.region, width=10000)
# Creates an InteractionSet object involving only bin pairs in the first and second regions.
# This provides a rapid subset of the functionality of squareCounts.
#
# written by Aaron Lun
# created 13 May 2017
# last modified 14 May 2017
{
    # Setting up the parameters
    width <- as.integer(width) 
    parsed <- .parseParam(param, bin=TRUE, width=width)
    chrs <- parsed$chrs
    frag.by.chr <- parsed$frag.by.chr
    discard <- parsed$discard
    cap <- parsed$cap
    bwidth <- parsed$bwidth
    bin.region <- parsed$bin.region
    bin.id <- parsed$bin.id
    bin.by.chr <- parsed$bin.by.chr

    # Checking the inputs.
    first.chr <- as.character(seqnames(first.region))
    second.chr <- as.character(seqnames(second.region))
    if (length(first.chr)!=1L) { stop("'first.region' should be of length 1") }
    if (length(second.chr)!=1L) { stop("'second.region' should be of length 1") }
    if (!(first.chr %in% chrs) || !(second.chr %in% chrs)) { 
        stop("anchor chromosome names not in cut site list") 
    }
    
    # Identifying the boxes that lie within our ranges of interest. We give it some leeway
    # to ensure that edges of the plot are retained.
    keep.frag.first <- keep.frag.second <- suppressWarnings(overlapsAny(bin.region, first.region)[bin.id])
    if (suppressWarnings(first.region!=second.region)) { 
        keep.frag.second <- suppressWarnings(overlapsAny(bin.region, second.region)[bin.id])
    }

    # Pulling out the read pair indices from each file, and checking whether chromosome names are flipped around.
    all.dex <- .loadIndices(file, chrs)
    flipped <- FALSE
    if (!is.null(all.dex[[first.chr]][[second.chr]])) {
        current <- .baseHiCParser(TRUE, file, first.chr, second.chr, 
            chr.limits=frag.by.chr, discard=discard, cap=cap, width=bwidth)[[1]]
        filter.a <- keep.frag.first
        filter.t <- keep.frag.second
        t.start <- bin.by.chr$first[[second.chr]]
        t.end <- bin.by.chr$last[[second.chr]]

    } else if (!is.null(all.dex[[second.chr]][[first.chr]])) { 
        current <- .baseHiCParser(TRUE, file, second.chr, first.chr, 
            chr.limits=frag.by.chr, discard=discard, cap=cap, width=bwidth)[[1]]
        filter.a <- keep.frag.second
        filter.t <- keep.frag.first
        t.start <- bin.by.chr$first[[first.chr]]
        t.end <- bin.by.chr$last[[first.chr]]
        flipped <- TRUE

    } else {
        # Just putting up some random values, if empty.
        filter.a <- filter.t <- logical(0) 
        t.start <- bin.by.chr$first[[second.chr]]
        t.end <- bin.by.chr$last[[second.chr]]
        current<-data.frame(anchor1.id=integer(0), anchor2.id=integer(0)) 
    }

    # Getting the read pairs around the area of interest.
    # Pick up reflection around diagonal (it's hard to conclusively define 
    # the anchor range acround the diagonal, so we just include everything).
    retain <- filter.a[current$anchor1.id] & filter.t[current$anchor2.id]
    if (first.chr==second.chr) { 
        retain <- retain | (filter.a[current$anchor2.id] & filter.t[current$anchor1.id]) 
    }

    # Collating them into counts and creating an output ISet.
    out <- .Call(cxx_count_patch, list(current[retain,]), bin.id, 1L, t.start, t.end)
    if (is.character(out)) { stop(out) }
    return(InteractionSet(list(counts=out[[3]]), metadata=List(param=param, width=width, flipped=flipped),
                          interactions=GInteractions(anchor1=out[[1]], anchor2=out[[2]], 
                                                     regions=bin.region, mode="reverse"))) 
}
