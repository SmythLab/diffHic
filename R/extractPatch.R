extractPatch <- function(file, param, first.region, second.region=first.region, width=10000, restrict.regions=FALSE)
# Creates an InteractionSet object involving only bin pairs in the first and second regions.
# This provides a rapid subset of the functionality of squareCounts.
#
# written by Aaron Lun
# created 13 May 2017
# last modified 18 November 2017
{
    # Setting up the parameters
    width <- as.integer(width) 
    first.chr <- as.character(seqnames(first.region))
    second.chr <- as.character(seqnames(second.region))
    if (length(first.chr)!=1L) { 
        stop("'first.region' should be of length 1") 
    }
    if (length(second.chr)!=1L) { 
        stop("'second.region' should be of length 1") 
    }

    # Constructing bins across the genome.
    param <- reform(param, restrict=cbind(first.chr, second.chr))
    is.dnase <- .isDNaseC(param)
    if (is.dnase) { 
        retainer <- c("anchor1.pos", "anchor2.pos", "anchor1.len", "anchor2.len")
        bin.out <- .createBins(param, width, restricted=restrict.regions)
    } else {
        retainer <- c("anchor1.id", "anchor2.id")
        bin.out <- .assignBins(param, width, restricted=restrict.regions)
    }
    bin.region <- bin.out$region
    bin.id <- bin.out$id
    bin.by.chr <- .splitByChr(bin.region)

    # Identifying the boxes that lie within our ranges of interest. 
    # Note: don't use overlapsAny(...)[bin.id], due to zeroes in bin.id when restrict.regions=TRUE.
    chrs <- seqlevelsInUse(bin.region)
    if (!(first.chr %in% chrs) || !(second.chr %in% chrs)) { 
        stop("anchor chromosome names not in cut site list") 
    }
    suppressWarnings(overlaps.first <- which(overlapsAny(bin.region, first.region))) 
    keep.frag.first <- keep.frag.second <- bin.id %in% overlaps.first
    if (suppressWarnings(first.region!=second.region)) { 
        suppressWarnings(overlaps.second <- which(overlapsAny(bin.region, second.region)))
        keep.frag.second <- bin.id %in% overlaps.second
    }

    # Pulling out the read pair indices from each file, and checking whether chromosome names are flipped around.
    loadfuns <- preloader(file, param=param, retain=retainer)
    flipped <- FALSE

    if (!is.null(loadfuns[[first.chr]][[second.chr]])) {
        current <- loadfuns[[first.chr]][[second.chr]][[1]]()
        t.start <- bin.by.chr$first[[second.chr]]
        t.end <- bin.by.chr$last[[second.chr]]
        if (is.dnase) { 
            a.start <- bin.by.chr$first[[first.chr]]
            a.end <- bin.by.chr$last[[first.chr]]
            current <- .binReads(current, width, a.start, t.start, a.end, t.end)
        }
        filter.a <- keep.frag.first
        filter.t <- keep.frag.second

    } else if (!is.null(loadfuns[[second.chr]][[first.chr]])) { 
        current <- loadfuns[[second.chr]][[first.chr]][[1]]()
        t.start <- bin.by.chr$first[[first.chr]]
        t.end <- bin.by.chr$last[[first.chr]]
        if (is.dnase) { 
            a.start <- bin.by.chr$first[[second.chr]]
            a.end <- bin.by.chr$last[[second.chr]]
            current <- .binReads(current, width, a.start, t.start, a.end, t.end)
        }
        filter.a <- keep.frag.second
        filter.t <- keep.frag.first
        flipped <- TRUE

    } else {
        # Just putting up some random values, if empty.
        filter.a <- filter.t <- logical(0) 
        t.start <- bin.by.chr$first[[second.chr]]
        t.end <- bin.by.chr$last[[second.chr]]
        current <- data.frame(anchor1.id=integer(0), anchor2.id=integer(0)) 
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
    return(InteractionSet(list(counts=out[[3]]), metadata=List(param=param, width=width, flipped=flipped),
                          interactions=GInteractions(anchor1=out[[1]], anchor2=out[[2]], 
                                                     regions=bin.region, mode="reverse"))) 
}
