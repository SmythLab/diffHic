###################################################################################################
# This tests the directionality index calculator.

Sys.setlocale(category="LC_COLLATE",locale="C")
chromos<-c(chrA=51, chrB=31)
source("simcounts.R")

dir.create("temp-domain")
dir1<-"temp-domain/1.h5"
dir2<-"temp-domain/2.h5"

suppressPackageStartupMessages(library(diffHic))
comp<-function(npairs1, npairs2, dist, cuts, restrict=NULL, cap=NA, span=5) {
    simgen(dir1, npairs1, chromos)
    simgen(dir2, npairs2, chromos)
    param <- pairParam(fragments=cuts, restrict=restrict, cap=cap)
    y <- squareCounts(c(dir1, dir2), param=param, width=dist, filter=1L)
    d <- domainDirections(c(dir1, dir2), param=param, width=dist, span=span)
    stopifnot(identical(regions(y), rowRanges(d)))

    collected.up <- collected.down <- matrix(0L, length(regions(y)), 2)
    for (chr in names(chromos)) {
        selected <- as.logical(seqnames(regions(y))==chr)

        for (lib in 1:2) {
            curmat <- as.matrix(inflate(y, chr, chr, sample=lib))
            up.counts <- down.counts <- integer(nrow(curmat))

            for (x in seq_len(nrow(curmat))) {
                indices <- x + seq_len(min(span, ncol(curmat)-x))
                up.counts[x] <- as.integer(sum(curmat[x, indices], na.rm=TRUE))
                indices <- x - seq_len(min(span, x-1))
                down.counts[x] <- as.integer(sum(curmat[indices, x], na.rm=TRUE))
            }

            collected.up[selected,lib] <- up.counts
            collected.down[selected,lib] <- down.counts
        }
    }

    ref.up <- assay(d, "up")
    ref.down <- assay(d, "down")
    dimnames(ref.up) <- dimnames(ref.down) <- NULL
    stopifnot(identical(collected.up, ref.up))
    stopifnot(identical(collected.down, ref.down))

    return(head(cbind(ref.up, ref.down)))
}

###################################################################################################

set.seed(100)
comp(200, 100, dist=10000, cuts=simcuts(chromos))
comp(200, 100, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(200, 100, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(200, 100, dist=10000, cuts=simcuts(chromos), span=2)
comp(200, 100, dist=10000, cuts=simcuts(chromos), span=10)
comp(200, 100, dist=5000, cuts=simcuts(chromos))
comp(200, 100, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(200, 100, dist=5000, cuts=simcuts(chromos), span=2)
comp(200, 100, dist=5000, cuts=simcuts(chromos), span=10)
comp(200, 100, dist=1000, cuts=simcuts(chromos))
comp(200, 100, dist=1000, cuts=simcuts(chromos), span=2)
comp(200, 100, dist=1000, cuts=simcuts(chromos), span=10)

comp(250, 200, dist=10000, cuts=simcuts(chromos))
comp(250, 200, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(250, 200, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(250, 200, dist=10000, cuts=simcuts(chromos), span=2)
comp(250, 200, dist=10000, cuts=simcuts(chromos), span=10)
comp(250, 200, dist=5000, cuts=simcuts(chromos))
comp(250, 200, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(250, 200, dist=5000, cuts=simcuts(chromos), span=2)
comp(250, 200, dist=5000, cuts=simcuts(chromos), span=10)
comp(250, 200, dist=1000, cuts=simcuts(chromos))
comp(250, 200, dist=1000, cuts=simcuts(chromos), span=2)
comp(250, 200, dist=1000, cuts=simcuts(chromos), span=10)

comp(500, 200, dist=10000, cuts=simcuts(chromos))
comp(500, 200, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(500, 200, dist=10000, cuts=simcuts(chromos, overlap=2))
comp(500, 200, dist=10000, cuts=simcuts(chromos), span=2)
comp(500, 200, dist=10000, cuts=simcuts(chromos), span=10)
comp(500, 200, dist=5000, cuts=simcuts(chromos))
comp(500, 200, dist=5000, cuts=simcuts(chromos, overlap=2))
comp(500, 200, dist=5000, cuts=simcuts(chromos), span=2)
comp(500, 200, dist=5000, cuts=simcuts(chromos), span=10)
comp(500, 200, dist=1000, cuts=simcuts(chromos))
comp(500, 200, dist=1000, cuts=simcuts(chromos), span=2)
comp(500, 200, dist=1000, cuts=simcuts(chromos), span=10)

##################################################################################################
# Cleaning up.

unlink("temp-domain", recursive=TRUE)

##################################################################################################
# End.


