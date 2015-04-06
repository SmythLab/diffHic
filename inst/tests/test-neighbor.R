###################################################################################################
# This tests the neighbor-counting code.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
suppressPackageStartupMessages(require(edgeR))

.getLimits <- function(x, flank, start, end) {
	lower.x <- x - flank
	upper.x <- x + flank
	if (lower.x < start) { upper.x <- upper.x + start - lower.x }
	if (upper.x > end) { lower.x <- lower.x + end - upper.x }
	lower.x <- max(start, lower.x)
	upper.x <- min(upper.x, end)
	return(c(lower.x, upper.x))
}

# Defining some odds and ends.

lower.left <- function(x) { 
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	out[nrow(x),1] <- FALSE
	out
}
upper.right <- function(x) { 
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	out[1, ncol(x)] <- FALSE
	out
}
upper.left <- function(x) { 
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	out[1,1] <- FALSE
	out	
}
lower.right <- function(x) { 
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	out[nrow(x), ncol(x)] <- FALSE
	out
}
all.but.middle <- function(x) {
	out <- matrix(TRUE, nrow=nrow(x), ncol=ncol(x))
	out[ceiling(length(out)/2)] <- FALSE
	out
}

comp <- function(npairs, chromos, flanking, prior=2) {
	flanking <- as.integer(flanking)

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
	bg <- enrichedPairs(data, flank=flanking, prior.count=prior)
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
				indices <- indices[keep(indices)]
				indices <- indices[indices > 0]
				indices <- indices[valid[indices]]

				# Computing the average across this quadrant.
				relevant.rows <- inter.space[indices]
				is.zero <- relevant.rows==0L			
				collected[quad] <- sum(rel.ab[relevant.rows[!is.zero]])/length(relevant.rows)
			}
#			print(sprintf("%i %i %.3f", ax, tx, collected[6]))
		
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

###################################################################################################

