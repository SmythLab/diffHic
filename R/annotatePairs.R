annotatePairs <- function(data.list, regions, rnames=names(regions), indices, ...) 
# Annotates clusters in 'indices' with 'rnames' based on overlapping 'regions' 
# for anchors of interactions in 'data.list'.
# 
# written by Aaron Lun
# created 13 January 2016    
{
    if (.notList(data.list)) { 
        if (missing(indices)) { 
            indices <- seq_along(data.list)
        }
        data.list <- list(data.list) 
    } else if (missing(indices)) { 
        stop("'indices' cannot be missing if 'data.list' is a list")
    }
    if (.notList(indices)) {
        indices <- list(indices) 
    }
    if (length(indices)!=length(data.list) || any(lengths(indices)!=lengths(data.list))) { 
        stop("length of elements in 'indices' and 'data.list' must be equal")
    }
    rnames <- as.character(rnames)
    if (length(regions)!=length(rnames)) {
        stop("length of 'regions' and 'rnames' must be the same")
    }

    # Identifying all overlapping features.
    ndata <- length(data.list)
    collected.anno1 <- collected.anno2 <- collected.index1 <- collected.index2 <- vector("list", ndata)
    for (x in seq_len(ndata)) { 
        curdata <- data.list[[x]]
        curdex <- indices[[x]]
        keep <- !is.na(curdex)
        curdata <- curdata[keep,]
        curdex <- curdex[keep]

        # Speeding up by overlapping just once.
        a1 <- anchors(curdata, type="first", id=TRUE)
        a2 <- anchors(curdata, type="second", id=TRUE)
        used <- union(a1, a2)
        olap <- findOverlaps(regions(curdata)[used], regions, ...)

        olap <- split(subjectHits(olap), queryHits(olap))
        all.olap <- rep(list(list()), length(regions(curdata)))
        all.olap[used[as.integer(names(olap))]] <- olap
        anno1 <- all.olap[a1] 
        anno2 <- all.olap[a2]

        collected.anno1[[x]] <- rnames[unlist(anno1)]
        collected.anno2[[x]] <- rnames[unlist(anno2)]
        collected.index1[[x]] <- curdex[rep(seq_along(anno1), lengths(anno1))]
        collected.index2[[x]] <- curdex[rep(seq_along(anno2), lengths(anno2))]
    }

    # Concatenating the strings.
    a1.anno <- sapply(split(unlist(collected.anno1), unlist(collected.index1)), .uniqConcat)
    a2.anno <- sapply(split(unlist(collected.anno2), unlist(collected.index2)), .uniqConcat)
    anchor1.anno <- anchor2.anno <- character(max(0L, sapply(indices, function(x) { if (length(x)) { return(max(x)) } else { return(0L) } })))
    anchor1.anno[as.integer(names(a1.anno))] <- a1.anno
    anchor2.anno[as.integer(names(a2.anno))] <- a2.anno
    return(list(anchor1=anchor1.anno, anchor2=anchor2.anno))
}

.uniqConcat <- function(x) { paste(unique(x), collapse=",") }

