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
		a<-anchors(y, id=TRUE, type="first")[x]
		t<-anchors(y, id=TRUE, type="second")[x]
		ref[a,]<-ref[a,]+ assay(y)[x,]
		ref[t,]<-ref[t,] + assay(y)[x,] 
	}

    out <- assay(frags)
    dimnames(out) <- NULL 
	if (!identical(ref, out)) { stop("mismatches in counts") }
	if (!identical(frags$totals, y$totals) || !identical(as.integer(colSums(assay(frags))), frags$totals*2L)) { 
		stop("mismatches in total counts") }
	if (!identical(regions(y), rowRanges(frags)))  { stop("mismatches in final regions") }

    if (!is.null(restrict)) { 
    	frag.alt <- marginCounts(c(dir1, dir2), param=param, width=dist, restrict.regions=TRUE)
        frag.sub <- frags[seqnames(rowRanges(frags)) %in% restrict,]
        if (!identical(assay(frag.sub), assay(frag.alt)) ||  !identical(rowRanges(frag.sub), rowRanges(frag.alt))) {
            stop("restrict.regions=TRUE doesn't work for marginCounts")
        }
    }

	return(head(assay(frags)))
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

