domainDirections <- function(files, param, width=50000, span=10)
# This function computes the directionality index for each bin in the genome,
# in order to find domain boundaries.
#
# written by Aaron Lun
# created 27 May 2015
# last modified 18 November 2017
{
    nlibs <- length(files)
    if (nlibs==0) { 
        stop("number of libraries must be positive")
    } else if (width < 0) { 
        stop("width must be a non-negative integer")
    } 
    width<-as.integer(width) 
    span <- as.integer(span)

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

    # Running through each pair of chromosomes.
    nlibs <- length(files)
    upcount <- downcount <- matrix(0L, length(bin.region), nlibs)

    loadfuns <- preloader(files, param=param, retain=retainer)
    for (chr in names(loadfuns)) {
        current <- loadfuns[[chr]]
        if (!(chr %in% names(current))) { 
            next 
        }
        curfuns <- current[[chr]]
        first.index <- bin.by.chr$first[[chr]]
        last.index <- bin.by.chr$last[[chr]]

        # Extracting counts.
        pairs <- vector("list", nlibs)
        for (lib in seq_len(nlibs)) { 
            cur.pairs <- curfuns[[lib]]()
            if (is.dnase) {
                cur.pairs <- .binReads(cur.pairs, width, first.index, first.index, last.index, last.index)
            }
            pairs[[lib]] <- cur.pairs
        }

        out <- .Call(cxx_directionality, pairs, bin.id, span, first.index, last.index)
        if (!length(out[[1]])) { 
            next 
        }
        pnts <- first.index:last.index
        downcount[pnts,] <- out[[1]]
        upcount[pnts,] <- out[[2]]
    }

    # Return an RSE with up and down counts.
    # No total counts, because we don't load every chromosome pair - might as well call totalCounts() externally if required.
    return(SummarizedExperiment(SimpleList(up=upcount, down=downcount), 
        bin.region, metadata=list(param=param, span=span, width=width)))                                      
}

