squareCounts <- function(files, param, width=50000, filter=1L, restrict.regions=FALSE)
# This function collates counts across multiple experiments to get the full set of results. It takes 
# a list of lists of lists of integer matrices (i.e. a list of the outputs of convertToInteractions) and
# then compiles the counts into a list object for output. 
#
# written by Aaron Lun
# some time ago
# last modified 17 November 2017
{
    nlibs <- length(files)
    if (nlibs==0L) {
        stop("number of libraries must be positive")
    } else if (width < 0) { 
        stop("width must be a non-negative integer")
    } 
    width <- as.integer(width) 
    filter <- as.integer(filter) 

    # Constructing bins across the genome.
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

    # Output vectors.
    full.sizes <- integer(nlibs)
    out.counts <- list(matrix(0L, 0, nlibs))
    out.a <- out.t <- list(integer(0))
    idex <- 1L

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

            # Switching to bin IDs for DNase-C data.
            if (is.dnase) {
                for (lib in seq_len(nlibs)) { 
                }
            }
            
            # Aggregating them in C++ to obtain count combinations for each bin pair.
            out <- .Call(cxx_count_patch, pairs, bin.id, filter, first.anchor2, last.anchor2)
            if (!length(out[[1]])) { 
                next 
            }

            # Storing counts and locations. 
            if (any(out[[1]] < out[[2]])) { stop("anchor1 ID should not be less than anchor2 ID") }
            out.a[[idex]] <- out[[1]]
             out.t[[idex]] <- out[[2]]
            out.counts[[idex]] <- out[[3]]
            idex<-idex+1L
        }
    }

    # Collating all the other results.
    out.a <- unlist(out.a)
    out.t <- unlist(out.t)
    out.counts <- do.call(rbind, out.counts)

    return(InteractionSet(list(counts=out.counts), colData=DataFrame(totals=full.sizes), 
        interactions=GInteractions(anchor1=out.a, anchor2=out.t, regions=bin.region, mode="reverse"), 
        metadata=List(param=param, width=width)))
}

## PROOF:
# Recall the enforcement of anchor1 >= anchor2. Bin pairs could technically be
# reflected around the diagonal, to ensure that all points are counted, e.g.,
# if a bin pair overlaps a point in (anchor2, anchor1) form but not in (anchor1,
# anchor2) form. However, this is not required, as shown below:
# 
# Consider a point (x, y) past the diagonal (i.e., no enforcement), where y >
# x. Assume that this point is covered by our bin pair. This implies that the
# anchor2 range of our bin pair `[te, ts]` includes 'y', and the anchor1 range
# `[ae, as]` includes 'x'. The anchor1 range must be above the anchor2 range
# (i.e., as >= ts, ae >= te). If ts <= y <= te and as <= x <= ae, then you can
# fiddle with this to obtain ts <= x <= te and as <= y <= ae (as x < y), i.e.,
# the point (y, x) is also covered. So, we're guaranteed to have already
# counted anything past the diagonal, meaning that reflection is not needed.

