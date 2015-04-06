# This script tests the capability of the filtering methods, on filtration 
# of read pairs during count loading, e.g. based on gap distances. We only
# check the .baseHiCParser function here, as the interface is common.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
suppressPackageStartupMessages(require(rhdf5))
source("simcounts.R")

dir.create("temp-filt")
dir1<-"temp-filt/1.h5"
dir2<-"temp-filt/2.h5"
dir1x<-"temp-filt/1b.h5"
dir2x<-"temp-filt/2b.h5"

filtsim <- function(npairs1, npairs2, chromos, overlap=4, min.ingap=NA, min.outgap=NA, max.frag=NA, discard.param=NULL, cap=NA) {
	simgen(dir1, npairs1, chromos)
	simgen(dir2, npairs2, chromos)
	cuts <- simcuts(chromos, overlap=overlap)
	augmentsim(dir1, cuts)
	augmentsim(dir2, cuts)
	cap <- as.integer(cap)

	# Pruning.
	totes1 <- prunePairs(dir1, pairParam(fragments=cuts), file.out=dir1x, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)
	totes2 <- prunePairs(dir2, pairParam(fragments=cuts), file.out=dir2x, min.inward=min.ingap, min.outward=min.outgap, max.frag=max.frag)

	# Setting up parameters.
	min.ingap <- as.integer(min.ingap)
	min.outgap <- as.integer(min.outgap)
	max.frag <- as.integer(max.frag)
	if (is.null(discard.param)) { 
		discard <- NULL
	} else {
		discard <- makeDiscard(discard.param[1], discard.param[2], seqlengths(cuts))
	}
	xdiscard <- diffHic:::.splitDiscards(discard)
	lost.frag <- lost.in <- lost.out <- lost.disc <- 0L

	combo <- c(dir1, dir2, dir1x, dir2x)
	stuff <- diffHic:::.loadIndices(combo, seqlevels(cuts))
	for (ax in names(stuff)) { 
		current <- stuff[[ax]]
		for (tx in names(current)) { 

			everything <- list()
			for (d in 1:2) {
				if (current[[tx]][d]) { 
					collected <- h5read(combo[d], file.path(ax, tx))
				} else {
					collected <- data.frame(anchor.id=integer(0),
						target.id=integer(0), anchor.pos=integer(0), target.pos=integer(0),
						anchor.len=integer(0), target.len=integer(0))
				}
				mod.d <- d+2
				if (current[[tx]][mod.d]) { 
					counter <- h5read(combo[mod.d], file.path(ax, tx))
				} else {
					counter <- data.frame(anchor.id=integer(0),
						target.id=integer(0), anchor.pos=integer(0), target.pos=integer(0),
						anchor.len=integer(0), target.len=integer(0))
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
				for (x in 1:ncol(counter)) { 
					attributes(mod.col[[x]]) <- attributes(counter[[x]]) <- NULL 
				}
				stopifnot(identical(counter, mod.col)) 
				
				# Checking discards.
				disc.keep <- !logical(nrow(collected))
				if (!is.null(discard)) { 
					nokill <- !overlapsAny(GRanges(ax, IRanges(collected$anchor.pos, collected$anchor.pos+abs(collected$anchor.len)-1L)), discard, type="within") & 
						!overlapsAny(GRanges(tx, IRanges(collected$target.pos, collected$target.pos+abs(collected$target.len)-1L)), discard, type="within")
					lost.disc <- lost.disc + sum(!nokill)
					disc.keep <- nokill 
				}
				collected <- collected[disc.keep,1:2]

				# Comparing it to the cap.
				if (!is.na(cap)) { 
					is.diff <- c(TRUE, diff(collected$anchor.id)!=0L | diff(collected$target.id)!=0L)
					uniq.ids <- cumsum(is.diff)
					by.ids <- split(1:nrow(collected), uniq.ids)
					to.keep <- logical(nrow(collected))
					for (x in by.ids) {
						if (length(x)>cap) { x <- x[1:cap] } 
						to.keep[x] <- TRUE
					}
					collected <- collected[to.keep,]
				} 

				dim(collected$anchor.id) <- dim(collected$target.id) <- NULL
				everything[[d]] <- collected
			}
									
			# Running through them with every possible check.
			comparator <- diffHic:::.baseHiCParser(current[[tx]][1:2], combo[1:2], ax, tx, discard=xdiscard, cap=cap)
			for (d in 1:2) {
				stopifnot(identical(everything[[d]], comparator[[d]]))
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
filtsim(100, 100, chromos, 4, min.ingap=100)
filtsim(100, 100, chromos, 4, min.ingap=1000)
filtsim(100, 100, chromos, 4, min.outgap=100)
filtsim(100, 100, chromos, 4, min.outgap=1000)
filtsim(100, 100, chromos, 4, max.frag=100)
filtsim(100, 100, chromos, 4, max.frag=1000)

filtsim(200, 50, chromos, 2)
filtsim(200, 50, chromos, 2, min.ingap=100, min.outgap=1000)
filtsim(200, 50, chromos, 2, min.ingap=1000, min.outgap=100)
filtsim(200, 50, chromos, 2, min.outgap=100, max.frag=100)
filtsim(200, 50, chromos, 2, min.outgap=1000, max.frag=1000)
filtsim(200, 50, chromos, 2, max.frag=100, min.ingap=100)
filtsim(200, 50, chromos, 2, max.frag=1000, min.ingap=1000)

filtsim(200, 200, chromos, 2, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, min.ingap=100, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, min.ingap=1000, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, min.outgap=100, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, min.outgap=1000, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, max.frag=100, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, max.frag=1000, discard.param=c(10, 100))

filtsim(200, 200, chromos, 2, cap=1, min.ingap=200)
filtsim(200, 200, chromos, 2, cap=1, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, min.ingap=100, min.outgap=1000, cap=2)
filtsim(200, 200, chromos, 2, cap=2, max.frag=100)
filtsim(200, 200, chromos, 2, cap=5, max.frag=1000, min.ingap=1000)
filtsim(200, 200, chromos, 2, cap=5, min.outgap=100)
filtsim(200, 200, chromos, 2, ca=10, min.ingap=100, discard.param=c(10, 100))
filtsim(200, 200, chromos, 2, ca=10, discard.param=c(10, 100))

#########################################################################

unlink("temp-filt", recursive=TRUE)

#########################################################################
