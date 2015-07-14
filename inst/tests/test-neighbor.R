###################################################################################################
# This tests the neighbor-counting code.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
suppressPackageStartupMessages(require(edgeR))

# Defining some odds and ends.

lower.left <- function(x, exclude=0) { 
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	out[nrow(x)+(-exclude):0,1:(1+exclude)] <- FALSE
	out
}

all.but.middle <- function(x, exclude=0) {
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	midrow <- ceiling(nrow(x)/2) + (-exclude):exclude
	midrow <- midrow[midrow > 0 & midrow <= nrow(x)]
	midcol <- ceiling(ncol(x)/2) + (-exclude):exclude
	midcol <- midcol[midcol > 0 & midcol <= ncol(x)]
	out[midrow, midcol] <- FALSE
	out
}

comp <- function(npairs, chromos, flanking, exclude=0, prior=2) {
	flanking <- as.integer(flanking)
	exclude <- as.integer(exclude)

	nlibs <- 4L
	lambda <- 5
	nbins <- sum(chromos)
	all.pairs <- rbind(t(combn(nbins, 2)), cbind(1:nbins, 1:nbins))
	aid <- pmax(all.pairs[,1], all.pairs[,2])
	tid <- pmin(all.pairs[,1], all.pairs[,2])
   	npairs <- min(npairs, nrow(all.pairs))

	# Setting up some data.
	counts <- do.call(cbind, lapply(1:nlibs, FUN=function(x) { as.integer(rpois(npairs, lambda) + 1) }) )
	chosen <- sample(nrow(all.pairs), npairs)
	indices <- unlist(sapply(chromos, FUN=function(x) { 1:x }), use.names=FALSE)
	data <- DIList(counts=counts, anchors=aid[chosen], targets=tid[chosen],
		totals=rep(1e6, nlibs), regions=GRanges(rep(names(chromos), chromos), IRanges(indices, indices)))
	data@regions$nfrags <- rep(1:3, length.out=nbins)
	
	# Computing the reference enrichment value.
	bg <- enrichedPairs(data, flank=flanking, prior.count=prior, exclude=exclude)
	final.ref <- numeric(length(bg))

	# Sorting them by chromosome pairs.
	all.chrs <- as.character(seqnames(regions(data)))
	chr.pair <- paste0(all.chrs[data@anchors], ".", all.chrs[data@targets])
	by.chr.pair <- split(1:npairs, chr.pair)
	first.id <- lapply(split(1:nbins, all.chrs), FUN=min)

	for (cpair in names(by.chr.pair)) { 
		cur.pairs <- by.chr.pair[[cpair]]
		two.chrs <- strsplit(cpair, "\\.")[[1]]
		current <- data[cur.pairs,]
		rel.ab <- 2^(aveLogCPM(counts(current), lib.size=current$totals, prior.count=0) 
			+ log2(mean(current$totals)/1e6))

		# Setting up the interaction space.
		a.dex <- anchors(current, id=TRUE) - first.id[[two.chrs[1]]] + 1L
		t.dex <- targets(current, id=TRUE) - first.id[[two.chrs[2]]] + 1L
		alen <- chromos[[two.chrs[1]]]
		tlen <- chromos[[two.chrs[2]]]
		inter.space <- matrix(0L, nrow=alen, ncol=tlen)
		inter.space[(t.dex-1)*alen + a.dex] <- 1:nrow(current) # column major.
		valid <- matrix(TRUE, nrow=alen, ncol=tlen)
		
		# Checking if we're working on the same chromosome.
		if (two.chrs[1]==two.chrs[2]) { 
			valid[upper.tri(valid)] <- FALSE 
			starting.dex <- 1L
		} else {
			starting.dex <- 2L
		}

		output <- numeric(nrow(current))
		for (pair in 1:nrow(current)) {
			total.num <- 4L
			collected <- numeric(total.num)
			collected.n <- numeric(total.num)
			ax <- a.dex[pair]
			tx <- t.dex[pair]

			for (quad in starting.dex:total.num) {
				if (quad==1L) {
					cur.a <- ax - flanking:0
					cur.t <- tx + 0:flanking
					keep <- lower.left 
				} else if (quad==2L) {
					cur.a <- ax + (-flanking):flanking
					cur.t <- tx
					keep <- all.but.middle
				} else if (quad==3L) {
					cur.a <- ax
					cur.t <- tx + (-flanking):flanking
					keep <- all.but.middle
				} else if (quad==4L) {
					cur.a <- ax + (-flanking):flanking
					cur.t <- tx + (-flanking):flanking
					keep <- all.but.middle
				}
	
				# Selecting the relevant entries for the chosen quadrant.
				indices <- outer(cur.a, cur.t, FUN=function(x, y) { 
					out <- (y-1)*alen + x
					out[x > alen | x < 1 | y > tlen | y < 1] <- -1
					return(out)
				})
				indices <- indices[keep(indices, exclude)]
				indices <- indices[indices > 0]
				indices <- indices[valid[indices]]

				# Computing the average across this quadrant.
				relevant.rows <- inter.space[indices]
				is.zero <- relevant.rows==0L			
				collected[quad] <- sum(rel.ab[relevant.rows[!is.zero]])/length(relevant.rows)
				collected.n[quad] <- length(relevant.rows)
			}

#			if (exclude) { # Troubleshooting.
#				print(c(aid[pair], tid[pair]))
#				print(collected)
#				print(collected.n)
#			}
		
			output[pair] <- log2((rel.ab[pair]+prior)/(max(collected, na.rm=TRUE)+prior))
		}
		final.ref[cur.pairs] <- output
	}

	if (any(abs(bg-final.ref) > (0.001+abs(bg))*1e-6)) { stop("mismatch in relative enrichment values") }
	return(head(bg))
}

###################################################################################################
# Simulating.

set.seed(3427675)
comp(10, c(chrA=10), 5)
comp(100, c(chrA=10, chrB=30, chrC=20), 5)
comp(100, c(chrA=10, chrC=20), 5)
comp(100, c(chrA=10, chrB=5, chrC=20), 5)
comp(100, c(chrA=20, chrB=5), 5)

comp(100, c(chrA=10, chrB=30, chrC=20), 10)
comp(100, c(chrA=10, chrC=20), 10)
comp(100, c(chrA=10, chrB=5, chrC=20), 10)
comp(100, c(chrA=20, chrB=10), 10)

comp(200, c(chrA=10, chrB=30, chrC=20), 3)
comp(200, c(chrA=10, chrC=20), 3)
comp(200, c(chrA=10, chrB=5, chrC=20), 3)
comp(200, c(chrA=20, chrB=3), 3)

comp(200, c(chrA=10, chrB=30, chrC=20), 1)
comp(200, c(chrA=10, chrC=20), 1)
comp(200, c(chrA=10, chrB=5, chrC=20), 1)
comp(200, c(chrA=20, chrB=5), 1)

comp(200, c(chrA=10, chrB=30, chrC=20), 3, exclude=1)
comp(200, c(chrA=10, chrC=20), 3, exclude=1)
comp(200, c(chrA=10, chrB=5, chrC=20), 3, exclude=1)
comp(200, c(chrA=20, chrB=5), 3, exclude=1)

###################################################################################################
# Same sort of simulation, but direct from read data, for neighbourCounts testing.

chromos<-c(chrA=51, chrB=31)
source("simcounts.R")

dir.create("temp-neighbor")
dir1<-"temp-neighbor/1.h5"
dir2<-"temp-neighbor/2.h5"

comp2 <- function(npairs1, npairs2, width, cuts, filter=1, flank=5, exclude=0, prior.count=2) {
	simgen(dir1, npairs1, chromos)
	simgen(dir2, npairs2, chromos)
	param <- pairParam(fragments=cuts)

	out <- neighbourCounts(c(dir1, dir2), param, width=width, filter=filter, flank=flank, prior.count=prior.count, 
		exclude=exclude)

	ref <- squareCounts(c(dir1, dir2), width=width, param, filter=1)
	keep <- rowSums(counts(ref)) >= filter
	enrichment <- enrichedPairs(ref, flank=flank, prior.count=prior.count, exclude=exclude)

	if (!identical(ref[keep,], out$interaction)) { stop("extracted counts don't match up") }
	if (any(abs(enrichment[keep] - out$enrichment) > 1e-6)) { stop("enrichment values don't match up") }
	return(head(enrichment[keep]))
}

set.seed(2384)
comp2(100, 50, 10000, cuts=simcuts(chromos))
comp2(100, 50, 10000, cuts=simcuts(chromos), filter=10)
comp2(100, 50, 10000, cuts=simcuts(chromos), flank=3)
comp2(100, 50, 10000, cuts=simcuts(chromos), prior.count=1)

comp2(50, 200, 5000, cuts=simcuts(chromos))
comp2(50, 200, 5000, cuts=simcuts(chromos), filter=10)
comp2(50, 200, 5000, cuts=simcuts(chromos), flank=3)
comp2(50, 200, 5000, cuts=simcuts(chromos), prior.count=1)

comp2(100, 200, 1000, cuts=simcuts(chromos))
comp2(100, 200, 1000, cuts=simcuts(chromos), filter=5)
comp2(100, 200, 1000, cuts=simcuts(chromos), flank=3)
comp2(100, 200, 1000, cuts=simcuts(chromos), prior.count=1)

comp2(10, 20, 1000, cuts=simcuts(chromos))
comp2(10, 20, 1000, cuts=simcuts(chromos), filter=5)
comp2(10, 20, 1000, cuts=simcuts(chromos), flank=3)
comp2(10, 20, 1000, cuts=simcuts(chromos), prior.count=1)

comp2(10, 20, 1000, cuts=simcuts(chromos), exclude=1)
comp2(50, 20, 1000, cuts=simcuts(chromos), exclude=1)
comp2(100, 50, 1000, cuts=simcuts(chromos), exclude=2)
comp2(50, 200, 1000, cuts=simcuts(chromos), exclude=2)

#####################################################################################################
# Cleaning up

unlink("temp-neighbor", recursive=TRUE)

#####################################################################################################
# End.


