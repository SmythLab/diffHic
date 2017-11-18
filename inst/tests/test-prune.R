# This script tests the capability of the filtering methods, on filtration 
# of read pairs during count loading, e.g. based on gap distances. 

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
suppressPackageStartupMessages(require(rhdf5))
source("simcounts.R")

dir.create("temp-filt")
dir1 <- "temp-filt/1.h5"
dir2 <- "temp-filt/2.h5"
dir1x <- "temp-filt/1a.h5"
dir2x <- "temp-filt/2a.h5"
dir1y <- "temp-filt/1b.h5"
dir2y <- "temp-filt/2b.h5"

filtsim <- function(npairs1, npairs2, chromos, overlap=4, min.ingap=NA, min.outgap=NA, max.frag=NA) { 
	simgen(dir1, npairs1, chromos)
	simgen(dir2, npairs2, chromos)
	cuts <- simcuts(chromos, overlap=overlap)
	augmentsim(dir1, cuts)
	augmentsim(dir2, cuts)

	# Straightforward pruning on min/max gap/frag.
    param <- pairParam(fragments=cuts)
  	totes1 <- prunePairs(dir1, param, file.out=dir1x, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)
	totes2 <- prunePairs(dir2, param, file.out=dir2x, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)

    # Setting up parameters.
	min.ingap <- as.integer(min.ingap)
	min.outgap <- as.integer(min.outgap)
	max.frag <- as.integer(max.frag)
	lost.frag <- lost.in <- lost.out <- lost.disc <- 0L

	combo <- c(dir1, dir2, dir1x, dir2x)
	stuff <- diffHic:::preloader(combo, param)
	for (ax in names(stuff)) { 
		current <- stuff[[ax]]
		for (tx in names(current)) { 

			everything <- list()
			for (i in 1:2) {
				if (!is.null(current[[tx]][[i]])) { 
					collected <- h5read(combo[i], file.path(ax, tx))
				} else {
					collected <- data.frame(anchor1.id=integer(0), anchor2.id=integer(0), 
                                            anchor1.pos=integer(0), anchor2.pos=integer(0),
                                            anchor1.len=integer(0), anchor2.len=integer(0))
				}
				mod.i <- i+2
				if (!is.null(current[[tx]][[mod.i]])) { 
					counter <- h5read(combo[mod.i], file.path(ax, tx))
				} else {
					counter <- data.frame(anchor1.id=integer(0), anchor2.id=integer(0), 
                                          anchor1.pos=integer(0), anchor2.pos=integer(0),
                                          anchor1.len=integer(0), anchor2.len=integer(0))
				}

				# Pruning the statistics.
				stats <- diffHic:::.getStats(collected, ax==tx, cuts)
				keep <- !logical(nrow(collected))
				if (!is.na(max.frag)) { 
					nokill <- stats$length <= max.frag
					lost.frag <- lost.frag + sum(!nokill)
					keep <- keep & nokill
				} 
				if (!is.na(min.outgap)) {
					nokill <- (ax!=tx | stats$orientation!=2L | stats$insert >= min.outgap)
					lost.out <- lost.out + sum(!nokill)
					keep <- keep & nokill 
				}
				if (!is.na(min.ingap)) { 
					nokill <- (ax!=tx | stats$orientation!=1L | stats$insert >= min.ingap)
					lost.in <- lost.in + sum(!nokill)
					keep <- keep & nokill 
				}
                mod.col <- collected[keep,]
				row.names(mod.col) <- NULL
				for (x in seq_len(ncol(counter))) { 
					attributes(mod.col[[x]]) <- attributes(counter[[x]]) <- NULL 
				}
				stopifnot(identical(counter, mod.col)) 
			}
		}
	}

    stopifnot(identical(lost.frag, totes1[["length"]]+totes2[["length"]]))	
	stopifnot(identical(lost.in, totes1[["inward"]]+totes2[["inward"]]))	
	stopifnot(identical(lost.out, totes1[["outward"]]+totes2[["outward"]]))	
	return(c(by.frag=lost.frag, by.in=lost.in, by.out=lost.out, by.disc=lost.disc))
}

#########################################################################

set.seed(38247)
chromos <- c(chrA=50, chrB=20, chrC=30)

filtsim(100, 100, chromos, 4)
filtsim(100, 100, chromos, 4, min.ingap=10000)
filtsim(100, 100, chromos, 4, min.ingap=50000)
filtsim(100, 100, chromos, 4, min.outgap=10000)
filtsim(100, 100, chromos, 4, min.outgap=50000)
filtsim(100, 100, chromos, 4, max.frag=100)
filtsim(100, 100, chromos, 4, max.frag=1000)

filtsim(200, 50, chromos, 2)
filtsim(200, 50, chromos, 2, min.ingap=10000, min.outgap=10000)
filtsim(200, 50, chromos, 2, min.ingap=50000, min.outgap=10000)
filtsim(200, 50, chromos, 2, min.outgap=10000, max.frag=100)
filtsim(200, 50, chromos, 2, min.outgap=10000, max.frag=1000)
filtsim(200, 50, chromos, 2, max.frag=100, min.ingap=10000)
filtsim(200, 50, chromos, 2, max.frag=1000, min.ingap=10000)

#########################################################################

addparamsim <- function(npairs1, npairs2, chromos, overlap=4, min.ingap=NA, min.outgap=NA, max.frag=NA, restrict=NULL, discard.param=NULL, cap=NA) { 
 	simgen(dir1, npairs1, chromos)
	simgen(dir2, npairs2, chromos)
	cuts <- simcuts(chromos, overlap=overlap)
	augmentsim(dir1, cuts)
	augmentsim(dir2, cuts)

    # Pruning with constant arguments.
    param <- pairParam(fragments=cuts)
    if (is.na(min.ingap) && is.na(min.outgap) && is.na(max.frag)) { 
        dir1x <- dir1
        dir2x <- dir2
    } else {
        prunePairs(dir1, param, file.out=dir1x, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)
        prunePairs(dir2, param, file.out=dir2x, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)
    }

   # Pruning with additional arguments.
    new.param <- pairParam(fragments=cuts, restrict=restrict, cap=cap)
	if (!is.null(discard.param)) { 
		discard <- makeDiscard(discard.param[1], discard.param[2], seqlengths(cuts))
        new.param <- reform(new.param, discard=discard)
	}
	totes1 <- prunePairs(dir1, new.param, file.out=dir1y, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)
	totes2 <- prunePairs(dir2, new.param, file.out=dir2y, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)

    # Checking that the counting behaves the same.
    ref <- squareCounts(c(dir1x, dir2x), new.param, width=1)
    out <- squareCounts(c(dir1y, dir2y), param, width=1)
    stopifnot(identical(interactions(ref), interactions(out)))
    stopifnot(identical(assay(ref), assay(out)))
    stopifnot(identical(ref$totals, out$totals))

    # Additional statistics check with getPairData on the pruned files.
    if (is.na(min.ingap) && is.na(min.outgap) && is.na(max.frag)) { 
        ref1 <- getPairData(dir1, new.param)
        out1 <- getPairData(dir1y, param)
        stopifnot(identical(ref1, out1))
        ref2 <- getPairData(dir2, new.param)
        out2 <- getPairData(dir2y, param)
        stopifnot(identical(ref2, out2))
    }

    return(totes1 + totes2)
}

set.seed(10893)
addparamsim(200, 200, chromos, 2) # for reference.

addparamsim(200, 200, chromos, 2, discard.param=c(10, 100))
addparamsim(200, 200, chromos, 2, discard.param=c(100, 200))
addparamsim(200, 200, chromos, 2, min.ingap=1000, discard.param=c(10, 100))
addparamsim(200, 200, chromos, 2, min.outgap=1000, discard.param=c(10, 100))
addparamsim(200, 200, chromos, 2, max.frag=100, discard.param=c(10, 100))

addparamsim(200, 200, chromos, 2, restrict="chrB")
addparamsim(200, 200, chromos, 2, restrict=cbind("chrA", "chrB"))
addparamsim(200, 200, chromos, 2, min.ingap=1000, restrict="chrB") 
addparamsim(200, 200, chromos, 2, min.outgap=1000, restrict="chrB") 
addparamsim(200, 200, chromos, 2, max.frag=100, restrict="chrB") 

addparamsim(200, 200, chromos, 2, cap=2)
addparamsim(200, 200, chromos, 2, cap=2, discard.param=c(10, 100))
addparamsim(200, 200, chromos, 2, cap=2, restrict="chrB")
addparamsim(200, 200, chromos, 2, min.ingap=1000, cap=2)
addparamsim(200, 200, chromos, 2, min.outgap=1000, cap=2)
addparamsim(200, 200, chromos, 2, max.frag=100, cap=2)

#########################################################################

unlink("temp-filt", recursive=TRUE)

#########################################################################
