###################################################################################################
# This tests the counting capabilities of various functions for DNase-C data.

Sys.setlocale(category="LC_COLLATE",locale="C")
dir.create("temp-dna")
file1 <- "temp-dna/1.h5"
file2 <- "temp-dna/2.h5"

suppressPackageStartupMessages(library(diffHic))

simDNA <- function(fout, chrs, npairs, rlen) { 
    r1 <- sample(length(chrs), npairs, replace=TRUE)
    r2 <- sample(length(chrs), npairs, replace=TRUE)
    p1 <- as.integer(runif(npairs, 1, chrs[r1] + 1))
    p2 <- as.integer(runif(npairs, 1, chrs[r2] + 1))
    l1 <- ifelse(rbinom(npairs, 1, 0.5)==1L, 1L, -1L)*rlen
    l2 <- ifelse(rbinom(npairs, 1, 0.5)==1L, 1L, -1L)*rlen

    savePairs(data.frame(anchor1.id=r1, anchor2.id=r2, anchor1.pos=p1, anchor2.pos=p2, anchor1.len=l1, anchor2.len=l2),
              fout, param=pairParam( GRanges(seqlengths=chrs) ))
    return(invisible(NULL))
}

comp <- function(chrs, npairs1, npairs2, dist, rlen=10, filter=1L, restrict=NULL, cap=NA) {
    simDNA(file1, chrs, npairs1, rlen)   
    simDNA(file2, chrs, npairs2, rlen)   

    # Output of squares.
    param <- pairParam( GRanges(seqlengths=chrs), restrict=restrict, cap=cap)
    y <- squareCounts(c(file1, file2), param=param, width=dist, filter=filter)

    # Reference. First, getting all bins.
    bin.coords <- list()
    for (i in names(chrs)) { 
        nbins <- ceiling(chrs[[i]]/dist)
        bin.ends <- pmin(chrs[[i]], seq_len(nbins)*dist)
        bin.starts <- c(1, head(bin.ends, -1)+1)
        bin.coords[[i]] <- GRanges(i, IRanges(bin.starts, bin.ends))
    }
    bin.offset <- c(0L, cumsum(lengths(bin.coords)))
    names(bin.coords) <- NULL
    suppressWarnings(bin.coords <- do.call(c, bin.coords))
    seqlengths(bin.coords) <- chrs
    bin.coords$nfrags <- 0L
    stopifnot(identical(bin.coords, regions(y)))

    # Now running through all bin pairs and assembling an InteractionSet.
    collected.isets <- list()
    collected.margins <- list()
    collected.totals <- list()
    for (f in 1:2) { 
        curf <- c(file1, file2)[f]
        fmat <- matrix(0L, length(bin.coords), length(bin.coords))
        total <- 0L

        for (i in seq_along(chrs)) { 
            for (j in seq_len(i)) {
                cur.i <- names(chrs)[i]
                cur.j <- names(chrs)[j]
                if (!is.null(restrict) && (!cur.i %in% restrict || !cur.j %in% restrict)) { next }
                curdat <- loadData(curf, cur.i, cur.j)
                total <- total+ nrow(curdat)
                
                p1 <- curdat$anchor1.pos + ifelse(curdat$anchor1.len > 0, 0L, -curdat$anchor1.len-1L)
                p1 <- pmin(p1, chrs[i])
                p2 <- curdat$anchor2.pos + ifelse(curdat$anchor2.len > 0, 0L, -curdat$anchor2.len-1L)
                p2 <- pmin(p2, chrs[j])
                b1 <- ceiling(p1/dist) + bin.offset[i]
                b2 <- ceiling(p2/dist) + bin.offset[j]

                for (x in seq_along(b1)) {
                    fmat[b1[x], b2[x]] <- fmat[b1[x], b2[x]] + 1L 
                    if (b1[x]!=b2[x]) { 
                        fmat[b2[x], b1[x]] <- fmat[b2[x], b1[x]] + 1L 
                    }
                }                
            }
        }

        cm <- ContactMatrix(fmat, seq_along(bin.coords), seq_along(bin.coords), regions=bin.coords)
        extractor <- upper.tri(fmat, diag=TRUE)
        is <- deflate(cm, extract=extractor)
        collected.isets[[f]] <- is
        collected.margins[[f]] <- as.integer(rowSums(fmat) + diag(fmat)) # diagonal gets counted twice.
        collected.totals[[f]] <- total
    }

    ref <- do.call(cbind, collected.isets)
    ref <- ref[rowSums(assay(ref)) >= filter,]
    interactions(ref) <- as(interactions(ref), "ReverseStrictGInteractions")
    storage.mode(assay(ref)) <- "integer"
    colnames(ref) <- NULL

    # Checking if interactions and counts are equal.
    m <- match(y, ref)
    stopifnot(identical(assay(ref)[m,], assay(y)))
    stopifnot(all(!is.na(m)))
    stopifnot(!anyDuplicated(m))
    stopifnot(nrow(y)==nrow(ref))
    
    # Checking the totals.
    stopifnot(identical(y$totals, as.integer(unlist(collected.totals))))
    if (is.null(restrict)) { 
        stopifnot(identical(y$totals, as.integer(c(npairs1, npairs2))))
    }
    totes <- totalCounts(c(file1, file2), param=param)
    stopifnot(identical(totes, y$totals))

    # Checking if the marginal counts... add up.
    mrg <- marginCounts(c(file1, file2), param=param, width=dist)
    out <- assay(mrg)
    dimnames(out) <- NULL
    stopifnot(identical(out, do.call(cbind, collected.margins)))
    stopifnot(identical(bin.coords, rowRanges(mrg)))

    # Checking that the neighborhood gives the same output.
    nbr <- neighborCounts(c(file1, file2), param=param, width=dist, filter=filter, flank=5)
    stopifnot(identical(assay(nbr), assay(y)))
    stopifnot(all(interactions(nbr)==interactions(y)))

    return(head(assay(y)))
}

set.seed(234711)

chrs <- c(chrA=1000, chrB=2000)
comp(chrs, 100, 200, 100, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 200, 200, 100, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 100, 200, 75, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 100, 200, 220, rlen=10, filter=1L, restrict=NULL, cap=NA)

comp(chrs, 500, 500, 75, rlen=10, filter=5L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 100, rlen=10, filter=5L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 150, rlen=10, filter=5L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 200, rlen=10, filter=5L, restrict=NULL, cap=NA)

comp(chrs, 500, 500, 75, rlen=50, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 100, rlen=50, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 150, rlen=50, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 200, rlen=50, filter=1L, restrict=NULL, cap=NA)

comp(chrs, 500, 500, 75, rlen=10, filter=1L, restrict="chrA", cap=NA)
comp(chrs, 500, 500, 100, rlen=10, filter=1L, restrict="chrA", cap=NA)
comp(chrs, 500, 500, 150, rlen=10, filter=1L, restrict="chrA", cap=NA)
comp(chrs, 500, 500, 200, rlen=10, filter=1L, restrict="chrA", cap=NA)

comp(chrs, 500, 500, 75, rlen=10, filter=1L, restrict=NULL, cap=5) # Should have no effect.
comp(chrs, 500, 500, 100, rlen=10, filter=1L, restrict=NULL, cap=5)
comp(chrs, 500, 500, 150, rlen=10, filter=1L, restrict=NULL, cap=5)
comp(chrs, 500, 500, 200, rlen=10, filter=1L, restrict=NULL, cap=5)

chrs <- c(chrA=1000, chrB=1000, chrC=1000) # Trying for more chromosomes.
comp(chrs, 500, 500, 75, rlen=10, filter=1L, restrict=NULL, cap=NA) 
comp(chrs, 500, 500, 100, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 150, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 200, rlen=10, filter=1L, restrict=NULL, cap=NA)

##################################################################################################
# Repeating the tests for connectCounts

simranges <- function(chrs, nranges, size.range=c(20, 100))
# Generates simulated ranges.
{
    ranges <- list()
    for (chr in names(chrs)) {
        chr.len <- chrs[[chr]]
        range.start <- round(runif(nranges, 1, chr.len))
        range.end <- pmin(chr.len, round(range.start + runif(nranges, size.range[1], size.range[2])))
        ranges[[chr]] <- GRanges(chr, IRanges(range.start, range.end))
    }               
    names(ranges) <- NULL
    suppressWarnings(ranges <- do.call(c, ranges))
    seqlengths(ranges) <- chrs
    return(ranges)  
}

concomp <- function(chrs, npairs1, npairs2, regions, rlen=10, filter=1L, restrict=NULL, cap=NA, seconds=NULL, type="any") {
    simDNA(file1, chrs, npairs1, rlen)   
    simDNA(file2, chrs, npairs2, rlen)   

    # Output of running connectCounts
    param <- pairParam( GRanges(seqlengths=chrs), restrict=restrict, cap=cap)
    y <- connectCounts(c(file1, file2), param=param, regions, filter=1L, second.regions=seconds, type=type)
    y.matches <- do.call(paste, anchors(y, id=TRUE))
    y.counts <- assay(y)

    if (filter>1L) {
        y.filt <- connectCounts(c(file1, file2), param=param, regions=regions, filter=filter, second.regions=seconds, type=type)
        stopifnot(isTRUE(all.equal(y.filt, y[rowSums(assay(y)) >= filter,])))
    }

    # Reference block - first, making the regions.
    regions$nfrags <- 0L
    regions$original <- seq_along(regions)
    if (!is.null(seconds)) {
        if (is.numeric(seconds)) {
            extras <- diffHic:::.getBinID(param$fragments, seconds)$region
        } else {
            extras <- seconds
            extras$nfrags <- 0L
        }
        regions$is.second <- FALSE
        extras$original <- seq_along(extras)
        extras$is.second <- TRUE
        suppressWarnings(regions <- c(regions, extras))
    }
    regions <- sort(regions)
    stopifnot(all(regions==regions(y)))
    stopifnot(identical(regions$is.second, regions(y)$is.second))
    stopifnot(identical(regions$original, regions(y)$original))

    # Setting up the combining function.
    secondary <- regions$is.second
    combFUN <- function(hits1, hits2) {
        ex1 <- rep(hits1, each=length(hits2))
        ex2 <- rep(hits2, length(hits1))
 
        if (!is.null(secondary)) {
            keep <- secondary[ex1]!=secondary[ex2]
            ex1 <- ex1[keep]
            ex2 <- ex2[keep]
        }

        out1 <- pmax(ex1, ex2)
        out2 <- pmin(ex1, ex2)
        unique(paste(out1, out2))
    }

    # Running through them and making sure they match up.
    for (f in 1:2) { 
        curf <- c(file1, file2)[f]
        collected <- list()

        for (i in seq_along(chrs)) { 
            for (j in seq_len(i)) {
                cur.i <- names(chrs)[i]
                cur.j <- names(chrs)[j]
                if (!is.null(restrict) && (!cur.i %in% restrict || !cur.j %in% restrict)) { next }
                curdat <- loadData(curf, cur.i, cur.j)
                
                r1 <- GRanges(cur.i, IRanges(curdat$anchor1.pos, width=abs(curdat$anchor1.len)))
                r2 <- GRanges(cur.j, IRanges(curdat$anchor2.pos, width=abs(curdat$anchor2.len)))
                olap1 <- findOverlaps(r1, regions, type=type)
                olap2 <- findOverlaps(r2, regions, type=type)

                for (x in seq_len(nrow(curdat))) {
                    colap1 <- subjectHits(olap1)[queryHits(olap1)==x]
                    colap2 <- subjectHits(olap2)[queryHits(olap2)==x]
                    if (!length(colap1) || !length(colap2)) next
                    collected[[length(collected)+1L]] <- combFUN(colap1, colap2)
                }
            }
        }
        
        collected <- unlist(collected)
        col.counts <- table(collected)
        m <- match(names(col.counts), y.matches)
        stopifnot(all(!is.na(m)))
        y.counts[m,f] <- y.counts[m,f] - col.counts
    }
    stopifnot(all(y.counts==0L))

    return(head(assay(y)))
}

set.seed(1012)

chrs <- c(chrA=1000, chrB=2000)
regs <- simranges(chrs, 10)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA)
concomp(chrs, 100, 100, regs, rlen=50L, filter=1L, restrict=NULL, cap=NA)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict="chrB", cap=NA)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=1L)
concomp(chrs, 100, 100, regs, rlen=10L, filter=2L, restrict=NULL, cap=NA)

# Checking other types of overlaps
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA, type="within")
concomp(chrs, 100, 100, regs, rlen=50L, filter=1L, restrict=NULL, cap=NA, type="within")
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict="chrB", cap=NA, type="within")
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=1L, type="within")
concomp(chrs, 100, 100, regs, rlen=10L, filter=2L, restrict=NULL, cap=NA, type="within")

# Checking with bigger ranges.
regs <- simranges(chrs, 10, size.range=c(100, 500))
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA)
concomp(chrs, 100, 100, regs, rlen=50L, filter=1L, restrict=NULL, cap=NA)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict="chrB", cap=NA)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=1L)
concomp(chrs, 100, 100, regs, rlen=10L, filter=2L, restrict=NULL, cap=NA)

# Checking with more ranges.
chrs <- c(chrA=1000, chrB=1000, chrC=1000)
regs <- simranges(chrs, 20) 
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA)
concomp(chrs, 100, 100, regs, rlen=50L, filter=1L, restrict=NULL, cap=NA)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict="chrB", cap=NA)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=1L)
concomp(chrs, 100, 100, regs, rlen=10L, filter=2L, restrict=NULL, cap=NA)

# Checking with secondary regions.
chrs <- c(chrA=1000, chrB=2000)
regs <- simranges(chrs, 10)
regs2 <- simranges(chrs, 10)

concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA, seconds=regs2)
concomp(chrs, 100, 100, regs, rlen=50L, filter=1L, restrict=NULL, cap=NA, seconds=regs2)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict="chrB", cap=NA, seconds=regs2)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=1L, seconds=regs2)
concomp(chrs, 100, 100, regs, rlen=10L, filter=2L, restrict=NULL, cap=NA, seconds=regs2)

concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA, seconds=100)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA, seconds=200)
concomp(chrs, 100, 100, regs, rlen=10L, filter=1L, restrict=NULL, cap=NA, seconds=500)

##################################################################################################
# Miscellaneous bits and pieces, to check they work properly with DNase-C data.

# boxPairs.

set.seed(34234)
chrs <- c(chrA=1000, chrB=2000)
simDNA(file1, chrs, 1000, 10)   
simDNA(file2, chrs, 1000, 10)  

param <- pairParam( GRanges(seqlengths=chrs))
y1 <- squareCounts(c(file1, file2), param=param, width=100, filter=1L)
y2 <- squareCounts(c(file1, file2), param=param, width=200, filter=1L)

checkFUN <- function(y, box) {
    olap <- findOverlaps(y, box, type="within")
    stopifnot(identical(queryHits(olap), seq_len(nrow(y))))
    return(subjectHits(olap))
}

out <- boxPairs(y1, y2)
out$interactions
stopifnot(identical(out$indices[[1]], checkFUN(y1, out$interactions)))
stopifnot(identical(out$indices[[2]], checkFUN(y2, out$interactions)))

out <- boxPairs(y1, y2, reference=400)
out$interactions
stopifnot(identical(out$indices[[1]], checkFUN(y1, out$interactions)))
stopifnot(identical(out$indices[[2]], checkFUN(y2, out$interactions)))

# getArea

head(getArea(y1))
head(getArea(y1, bp=FALSE))

# Plotting.

pdf("temp-dna/out.pdf")
plotPlaid(file1, param, first.region=GRanges("chrA", IRanges(1, 100)), 
          second.region=GRanges("chrA", IRanges(1, 200)), width=50, diag=TRUE)
plotPlaid(file1, param, first.region=GRanges("chrA", IRanges(1, 100)), 
          second.region=GRanges("chrB", IRanges(1, 200)), width=50, diag=TRUE)
rotPlaid(file1, param, region=GRanges("chrA", IRanges(1, 200)), width=50)
rotPlaid(file1, param, region=GRanges("chrB", IRanges(1, 200)), width=50)
dev.off()

##################################################################################################
# Cleaning up.

unlink("temp-dna", recursive=TRUE)

##################################################################################################
# End.
