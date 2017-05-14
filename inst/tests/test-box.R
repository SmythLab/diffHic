###################################################################################################
# This tests the bin square summarization method in diffHic.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
source("simcounts.R")

chromos<-c(chrA=50, chrB=100, chrC=80)
comp <- function(reference, widths, minbox=FALSE) {
	cutted <- simcuts(chromos, min=20, max=200, overlap=4)
	collected <- list()
	for (w in widths) {
		output <- list()
		bindata <- diffHic:::.assignBins(cutted, w)
		n <- round(runif(1, 20, 100))
		all.a <- as.integer(runif(n, 1, length(bindata$region)))
		all.t <- as.integer(runif(n, 1, all.a))

		oname <- paste0("w", w)
		collected[[oname]] <- InteractionSet(list(counts=matrix(0L, nrow=n, ncol=1)), 
            colData=DataFrame(totals=100), 
            GInteractions(anchor1=all.a, anchor2=all.t, regions=bindata$region, mode="reverse"), 
			metadata=List(param=pairParam(fragments=cutted)))
	}	

	output<- do.call(boxPairs, c(collected, reference=reference, minbox=minbox))
	stopifnot(length(output$indices)==length(widths))
	for (x in 1:length(output$indices)) { 
		curdex <- output$indices[[x]]
		curlist <- collected[[x]]

		# Checking that each bin pair is truly nested within its reported parent.
		parent.a <- anchors(output$interactions[curdex], type="first")
		parent.t <- anchors(output$interactions[curdex], type="second")
		current.a <- anchors(curlist, type="first")
		current.t <- anchors(curlist, type="second")
		
		if (! all(seqnames(parent.a)==seqnames(current.a) & start(parent.a) <= start(current.a) & end(parent.a) >= end(current.a)) ) { stop("anchor ranges not nested in parent") }
		if (! all(seqnames(parent.t)==seqnames(current.t) & start(parent.t) <= start(current.t) & end(parent.t) >= end(current.t)) ) { stop("target ranges not nested in parent") }
	}
	
	if (minbox) { 
		# Checking that the minimum bounding box was correctly assigned.
		all.anchor1 <- all.anchor2 <- list()
		for (x in 1:length(output$indices)) { 
			all.anchor1[[x]] <- anchors(collected[[x]], type="first")
			all.anchor2[[x]] <- anchors(collected[[x]], type="second")
		}
		gathered.a <- unlist(range(split(do.call(c, all.anchor1), unlist(output$indices))))
		gathered.t <- unlist(range(split(do.call(c, all.anchor2), unlist(output$indices))))
		names(gathered.a) <- names(gathered.t) <- NULL
		if (!identical(gathered.a, anchors(output$interactions, type="first")) || 
            !identical(gathered.t, anchors(output$interactions, type="second"))) { stop("mismatch in bounding boxes") }
	}	

    nocom <- do.call(boxPairs, c(collected, reference=reference, minbox=minbox, index.only=TRUE))
    if (!identical(nocom, output$indices)) { stop("results aren't the same with index.only=TRUE") }
	return(lapply(output$indices, FUN=head))
}

set.seed(74653812)
comp(1000, c(100))
comp(1000, c(100, 500))
comp(1000, c(500))
comp(1000, c(500, 1000))
comp(1000, c(100, 500, 1000))

comp(150, c(10))
comp(150, c(10, 50))
comp(150, c(50))
comp(150, c(50, 150))
comp(150, c(10, 50, 150))

comp(500, c(50))
comp(500, c(250))
comp(500, c(50, 250))

comp(minbox=TRUE, 1000, c(100))
comp(minbox=TRUE, 1000, c(100, 500))
comp(minbox=TRUE, 1000, c(500))
comp(minbox=TRUE, 1000, c(500, 1000))
comp(minbox=TRUE, 1000, c(100, 500, 1000))

comp(minbox=TRUE, 150, c(10))
comp(minbox=TRUE, 150, c(10, 50))
comp(minbox=TRUE, 150, c(50))
comp(minbox=TRUE, 150, c(50, 150))
comp(minbox=TRUE, 150, c(10, 50, 150))

comp(minbox=TRUE, 500, c(50))
comp(minbox=TRUE, 500, c(250))
comp(minbox=TRUE, 500, c(50, 250))

####################################################################################################
# End.
