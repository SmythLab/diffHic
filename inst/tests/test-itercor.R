####################################################################################################
# Tests the iterative correction script.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
suppressPackageStartupMessages(require(edgeR))
	
comp<- function(npairs, nfrags, nlibs, lambda=5, dispersion=0.05, winsorize=0.02, discard=0.02, locality=1) {
	all.pairs <- rbind(t(combn(nfrags, 2)), cbind(1:nfrags, 1:nfrags))
	all.pairs <- data.frame(anchor.id=all.pairs[,2], target.id=all.pairs[,1])	
	npairs <- min(npairs, nrow(all.pairs))
	counts <- do.call(cbind, lapply(1:nlibs, FUN=function(x) { rpois(npairs, lambda) }) )
	chosen <- sample(nrow(all.pairs), npairs)
	data <- DIList(counts=counts, anchors=all.pairs$anchor.id[chosen], 
		targets=all.pairs$target.id[chosen], totals=rep(1, nlibs), 
		regions=GRanges(sort(sample(c("chrA", "chrB", "chrC"), nfrags, replace=TRUE)),
			IRanges(1:nfrags, 1:nfrags)))
	
	# Constructing the values.	
	actual.mat<-matrix(0, nfrags, nfrags)
	is.filled <- matrix(FALSE, nfrags, nfrags)
	ave.count <- exp(mglmOneGroup(counts, offset=numeric(nlibs), dispersion=dispersion))
	for (x in 1:nrow(data)) { 
		if (ave.count[x] < 1e-6) { next } # As zeros get removed.
		a<-data@anchors[x]
		t<-data@targets[x]
		actual.mat[a,t] <- ave.count[x]
		is.filled[a,t] <- TRUE
		if (a!=t) { 
			actual.mat[t,a] <- ave.count[x] 
			is.filled[t,a] <- TRUE
		}
	}

	# Negating local interations.
	if (locality >= 0L){
 		per.chr <- split(1:nfrags, as.integer(seqnames(regions(data))))
		for (curloc in 0:locality) {
			failed <- 0L
			for (chr in per.chr) {
				if (length(chr)<=curloc) { 
					failed <- failed + 1L
					next 
				}
				current <- chr[1:(length(chr) - curloc)]
				deleted <- (current - 1) * nfrags + current + curloc
				actual.mat[deleted] <- 0
				is.filled[deleted] <- FALSE
				
				deleted <- (current + curloc - 1) * nfrags + current
				actual.mat[deleted] <- 0
				is.filled[deleted] <- FALSE
			}
			if (failed==length(per.chr)) { break }
		}		
	}
	
	# Winsorizing.
	temp.mat <- actual.mat[lower.tri(actual.mat, diag=TRUE)]
	is.nonzero <- temp.mat>1e-6
	winsor.val <- max(temp.mat[is.nonzero][sum(is.nonzero) - rank(temp.mat[is.nonzero], ties="first") + 1L > sum(is.nonzero)*winsorize])
	actual.mat[actual.mat > winsor.val] <- winsor.val

	# Running the reference, and checking that the right number of low fragments are removed.
	# Can't do it directly, because sorting might not be consistent between R and C++.
	iters <- 50
	test <- correctedContact(data, dispersion=dispersion, winsor=winsorize, ignore=discard, 
			iterations=iters, exclude.local=locality)
	to.discard <- is.na(test$bias)
	frag.sum <- rowSums(actual.mat) 
	nempty.frags <- rowSums(is.filled) > 0L
	gap <- sum(to.discard & nempty.frags) - sum(nempty.frags) * discard
	if (gap >= 1e-6 || gap < -1 - 1e-6) { stop("number discarded is not consistent with that requested") }
	if (any(to.discard)) { stopifnot(max(frag.sum[to.discard]) <= min(frag.sum[!to.discard]) + 1e-6) }

	# Discarding those that have a low sum of counts.
	actual.mat[to.discard,] <- 0
	actual.mat[,to.discard] <- 0

	# Iterative correction.
	bias<-rep(1, nfrags)
	for (i in 1:iters) {
		additional<-sqrt(rowSums(actual.mat))
		bias <- bias*additional
		additional[additional==0]<-1
		actual.mat<-t(t(actual.mat/additional)/additional)
	}
	
	# Comparing to the reference implementation. We use a fairly gentle threshold for differences,
	# due to the iterative nature of things (and numerical instability and so forth).
	is.okay <- !to.discard
	if (any(abs(test$bias[is.okay]-bias[is.okay]) > 1e-6 * bias[is.okay])) { stop("biases do not match up") }
	return(head(bias))
}

####################################################################################################

set.seed(0)

# Varying the number of fragments.

comp(100, 20, 2, discard=0.1)
comp(100, 30, 2, discard=0.1)
comp(100, 40, 2, discard=0.1)

comp(100, 20, 2, winsor=0.05)
comp(100, 30, 2, winsor=0.1)
comp(100, 40, 2, winsor=0.01)

comp(100, 20, 2, locality=0)
comp(100, 30, 2, locality=2)
comp(100, 40, 2, locality=1)

# Trying with fewer reads.

#debug(comp)
comp(10, 20, 2, discard=0.1)
comp(10, 30, 2, discard=0.1)
comp(10, 40, 2, discard=0.1)

comp(20, 20, 2, winsor=0.05)
comp(20, 30, 2, winsor=0.1)
comp(20, 40, 2, winsor=0.01)

comp(10, 20, 2, locality=0)
comp(10, 30, 2, locality=2)
comp(10, 40, 2, locality=1)

# Trying with fewer libraries.

comp(50, 20, 1, discard=0.1)
comp(50, 30, 1, discard=0.1)
comp(50, 40, 1, discard=0.1)

comp(50, 20, 1, winsor=0.05)
comp(50, 30, 1, winsor=0.1)
comp(50, 40, 1, winsor=0.01)

comp(50, 20, 1, locality=0)
comp(50, 30, 1, locality=2)
comp(50, 40, 1, locality=1)

# Trying with no special attention.
comp(50, 20, 1, discard=0, winsor=0, locality=-1)
comp(50, 50, 1, discard=0, winsor=0, locality=-1)
comp(50, 20, 2, discard=0, winsor=0, locality=-1)
comp(50, 20, 2, discard=0, winsor=0, locality=1000)
comp(100, 20, 2, discard=0, winsor=0, locality=1000)

####################################################################################################
# End.
