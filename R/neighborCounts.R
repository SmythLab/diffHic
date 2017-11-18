neighborCounts <- function(files, param, width=50000, filter=1L, flank=NULL, exclude=NULL)
# This does the same thing as squareCounts, except that it simultaneously computes the 
# filter statistic for each extracted bin pair. This has lower memory requirements as
# it doesn't need to hold the entire `filter=1` output in memory at once.
#
# written by Aaron Lun
# created 21 May 2015
# last modified 18 November 2017
{
    nlibs <- length(files)
    if (nlibs==0L) {
        stop("number of libraries must be positive")
    } else if (width < 0) { 
        stop("width must be a non-negative integer")
    } 
    width<-as.integer(width) 
    filter<-as.integer(filter) 

    # Constructing bins across the genome.
    is.dnase <- .isDNaseC(param)
    if (is.dnase) { 
        retainer <- c("anchor1.pos", "anchor2.pos", "anchor1.len", "anchor2.len")
        bin.out <- .createBins(param, width, restricted=FALSE)
    } else {
        retainer <- c("anchor1.id", "anchor2.id")
        bin.out <- .assignBins(param, width, restricted=FALSE)
    }
    bin.region <- bin.out$region
    bin.id <- bin.out$id
    bin.by.chr <- .splitByChr(bin.region)

    # Output vectors.
    out.counts <- list(matrix(0L, 0, nlibs))
    out.a <- out.t <- list(integer(0))
    full.sizes <- integer(nlibs)

    modes <- .neighbor_locales()
    neighbor.counts <- lapply(modes, FUN=function(x) list(matrix(0L, 0, nlibs)))
    neighbor.N <- lapply(modes, FUN=function(x) list(double(0)))
    idex <- 1L
    
    # Other stuff related to calculation of the neighborhood regions.    
    if (is.null(flank)) { 
        flank <- formals(enrichedPairs)$flank 
    }
    if (is.null(exclude)) { 
        exclude <- formals(enrichedPairs)$exclude 
    }
    flank <- as.integer(flank)
    exclude <- as.integer(exclude)
    if (flank <= 0L) { 
        stop("flank width must be a positive integer") 
    }
    if (exclude < 0L) { 
        stop("exclude width must be a non-negative integer") 
    }
    if (flank <= exclude) { 
        stop("exclude width must be less than the flank width") 
    }

    # Running through each pair of chromosomes.
    loadfuns <- preloader(files, param=param, retain=retainer)
    for (anchor1 in names(loadfuns)) {
        current <- loadfuns[[anchor1]]
        first.anchor1 <- bin.by.chr$first[[anchor1]]
        last.anchor1 <- bin.by.chr$last[[anchor1]]

        for (anchor2 in names(current)) {
            curfuns <- current[[anchor2]] 
            first.anchor2 <- bin.by.chr$first[[anchor2]]
            last.anchor2 <- bin.by.chr$last[[anchor2]]

            # Extracting counts.
            pairs <- vector("list", nlibs)
            for (lib in seq_len(nlibs)) { 
                cur.pairs <- curfuns[[lib]]()
                if (is.dnase) {
                    cur.pairs <- .binReads(cur.pairs, width, first.anchor1, first.anchor2,
                                           last.anchor1, last.anchor2)
                }
                pairs[[lib]] <- cur.pairs
            }
            full.sizes <- full.sizes + sapply(pairs, FUN=nrow)

            # Aggregating counts in C++ to obtain count combinations for each bin pair.
            out <- .Call(cxx_count_background, pairs, bin.id, flank, exclude, filter, 
                first.anchor2, last.anchor2, first.anchor1, last.anchor1)
            if (!length(out[[1]])) { next }

            # Storing counts and locations. 
            if (any(out[[1]] < out[[2]])) { stop("anchor2 ID should not be less than anchor1 ID") }
            out.a[[idex]] <- out[[1]]
            out.t[[idex]] <- out[[2]]
            out.counts[[idex]] <- out[[3]]

            # Adding the neighbourhood statistics.
            for (m in seq_along(modes)) {
                neighbor.counts[[m]][[idex]] <- out[[4]][[m]]
                neighbor.N[[m]][[idex]] <- out[[5]][[m]]
            }
            idex <- idex+1L
        }
    }

    # Collating all the other results.
    out.a <- unlist(out.a)
    out.t <- unlist(out.t)
    out.counts <- do.call(rbind, out.counts)

    all.assays <- list(counts=out.counts)
    for (m in seq_along(modes)) { 
        all.assays[[modes[m]]] <- do.call(rbind, neighbor.counts[[m]])
    }

    out.IS <- InteractionSet(all.assays, colData=DataFrame(totals=full.sizes), 
        interactions=GInteractions(anchor1=out.a, anchor2=out.t, regions=bin.region, mode="reverse"), 
        metadata=List(param=param, width=width))
   
    n.names <- .neighbor_numbers() 
    for (m in seq_along(modes)) { 
        mcols(out.IS)[[n.names[m]]] <- unlist(neighbor.N[[m]])
    }
    return(out.IS)
}
