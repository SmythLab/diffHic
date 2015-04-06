###################################################################################################
# This tests the interaction counting capabilities of the marginal counter.

chromos<-c(chrA=51, chrB=31)
source("simcounts.R")

# We set up the comparison function to check our results. 

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
dir.create("temp-marg")
dir1<-"temp-marg/1.h5"
dir2<-"temp-marg/2.h5"

comp<-function(n1, n2, dist, cuts, restrict=NULL) {
	simgen(dir1, n1, chromos)
	simgen(dir2, n2, chromos)
	param <- pairParam(fragments=cuts, restrict=restrict)
	y<-squareCounts(c(dir1, dir2), param=param, width=dist, filter=1L)
	frags<-marginCounts(c(dir1, dir2), param=param, width=dist)
  
	n <- length(regions(y))
	ref <- matrix(0L, n, 2)
	for (x in 1:nrow(y)) {
		a<-y@anchors[x]
		t<-y@targets[x]
		ref[a,]<-ref[a,]+ counts(y)[x,]
		if (a!=t) { ref[t,]<-ref[t,]+counts(y)[x,] }
	}
	
	keep<-which(rowSums(ref)>0.5)
	if (!identical(ref[keep,], counts(frags))) { stop("mismatches in counts") }
	if (!identical(frags$totals, y$totals)) { stop("mismatches in total counts") }
	if (!identical(keep, frags@anchors) || !identical(keep, frags@targets)) { stop("mismatches in the regions to keep") }
	if (!identical(regions(y), regions(frags)))  { stop("mismatches in final regions") }
	return(head(counts(frags)))
}

###################################################################################################
# Checking a vanilla count.

set.seed(126857)
comp(20, 10, dist=10000, cuts=simcuts(chromos))
comp(20, 10, dist=10000, cuts=simcuts(chromos))
comp(20, 10, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(20, 10, dist=5000, cuts=simcuts(chromos))
comp(20, 10, dist=5000, cuts=simcuts(chromos, overlap=4))

# Repeating a couple of times.
comp(10, 10, dist=10000, cuts=simcuts(chromos))
comp(10, 10, dist=10000, cuts=simcuts(chromos))
comp(10, 10, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(10, 10, dist=5000, cuts=simcuts(chromos))
comp(10, 10, dist=5000, cuts=simcuts(chromos, overlap=4))

comp(10, 20, dist=10000, cuts=simcuts(chromos))
comp(10, 20, dist=10000, cuts=simcuts(chromos))
comp(10, 20, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(10, 20, dist=5000, cuts=simcuts(chromos))
comp(10, 20, dist=5000, cuts=simcuts(chromos, overlap=4))

###################################################################################################
# Another example, a bit more extreme with more overlaps.

comp(50, 20, dist=10000, cuts=simcuts(chromos))
comp(50, 20, dist=10000, cuts=simcuts(chromos))
comp(50, 20, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(50, 20, dist=5000, cuts=simcuts(chromos))
comp(50, 20, dist=5000, cuts=simcuts(chromos, overlap=4))

comp(30, 30, dist=10000, cuts=simcuts(chromos))
comp(30, 30, dist=10000, cuts=simcuts(chromos))
comp(30, 30, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(30, 30, dist=5000, cuts=simcuts(chromos))
comp(30, 30, dist=5000, cuts=simcuts(chromos, overlap=4))

###################################################################################################
# Another example which is the pinnacle of extremity.

comp(200, 100, dist=10000, cuts=simcuts(chromos))
comp(200, 100, dist=10000, cuts=simcuts(chromos))
comp(200, 100, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(200, 100, dist=5000, cuts=simcuts(chromos))
comp(200, 100, dist=5000, cuts=simcuts(chromos, overlap=4))

comp(50, 200, dist=10000, cuts=simcuts(chromos))
comp(50, 200, dist=10000, cuts=simcuts(chromos))
comp(50, 200, dist=10000, cuts=simcuts(chromos, overlap=4))
comp(50, 200, dist=5000, cuts=simcuts(chromos))
comp(50, 200, dist=5000, cuts=simcuts(chromos, overlap=4))

###################################################################################################
# Adding some restriction.

comp(20, 10, dist=10000, cuts=simcuts(chromos), restrict="chrA")
comp(20, 10, dist=10000, cuts=simcuts(chromos), restrict="chrB")
comp(20, 10, dist=10000, cuts=simcuts(chromos, overlap=4), restrict="chrA")
comp(20, 10, dist=5000, cuts=simcuts(chromos), restrict="chrB")

##################################################################################################
# Cleaning up.

unlink("temp-marg", recursive=TRUE)

##################################################################################################
# End.

