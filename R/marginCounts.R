marginCounts <- function(files, param, width=50000)
# Gets the marginal counts i.e. sum of counts for each bin or region.
# This is useful to determine the `genomic coverage' of each region,
# based on the number of Hi-C read pairs involving that region.
#
# written by Aaron Lun
# some time ago.
# last modified 18 November 2017
{
    nlibs <- length(files)
    width <- as.integer(width)
    if (width < 0) { stop("width must be a non-negative integer") }

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

    # Setting up output elements.
    total.counts <- matrix(0L, length(bin.region), nlibs)
    full.sizes <- integer(nlibs)

    # Running through each pair of chromosomes.
    loadfuns <- preloader(files, param=param, retain=retainer)
    for (anchor1 in names(loadfuns)) {
        current <- loadfuns[[anchor1]]
        first.anchor1 <- bin.by.chr$first[[anchor1]]
        last.anchor1 <- bin.by.chr$last[[anchor1]]
        nabins <- last.anchor1 - first.anchor1 + 1L
        keep.a <- first.anchor1:last.anchor1

        for (anchor2 in names(current)) {
            curfuns <- current[[anchor2]]
            first.anchor2 <- bin.by.chr$first[[anchor2]]
            last.anchor2 <- bin.by.chr$last[[anchor2]]
            ntbins <- last.anchor2 - first.anchor2 + 1L
            keep.t <- first.anchor2:last.anchor2

            # Aggregating them for each library.
            for (lib in seq_len(nlibs)) {
                cur.pairs <- curfuns[[lib]]()
                
                # Switching to bin IDs for DNase-C data.
                if (is.dnase) {
                    cur.pairs <- .binReads(cur.pairs, width, first.anchor1, first.anchor2,
                                           last.anchor1, last.anchor2)
                }

                a.counts <- tabulate(bin.id[cur.pairs$anchor1.id]-first.anchor1+1L, nbins=nabins)
                total.counts[keep.a,lib] <- total.counts[keep.a,lib] + a.counts
                t.counts <- tabulate(bin.id[cur.pairs$anchor2.id]-first.anchor2+1L, nbins=ntbins)
                total.counts[keep.t,lib] <- total.counts[keep.t,lib] + t.counts
                full.sizes[lib] <- full.sizes[lib] + nrow(cur.pairs)
            }
        }    
    }
    
    # Aggregating all elements.
    return(SummarizedExperiment(list(counts=total.counts), colData=DataFrame(totals=full.sizes), 
            rowRanges=bin.region, metadata=List(param=param)))
}

