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
		bindata <- diffHic:::.getBinID(cutted, w)
		n <- round(runif(1, 20, 100))
		all.a <- as.integer(runif(n, 1, length(bindata$region)))
		all.t <- as.integer(runif(n, 1, all.a))

		oname <- paste0("w", w)
		collected[[oname]] <- DIList(counts=matrix(0, nrow=n, ncol=1), totals=0, 
			anchors=all.a, targets=all.t, regions=bindata$region, 
			exptData=List(param=pairParam(fragments=cutted)))
	}	

	output<- do.call(boxPairs, c(collected, reference=reference, minbox=minbox))
	stopifnot(length(output$indices)==length(widths))
	for (x in 1:length(output$indices)) { 
		curdex <- output$indices[[x]]
		curlist <- collected[[x]]

		# Checking that each bin pair is truly nested within its reported parent.
		parent.a <- output$anchors[curdex]
		parent.t <- output$targets[curdex]
		current.a <- anchors(curlist)
		current.t <- targets(curlist)
		
		if (! all(seqnames(parent.a)==seqnames(current.a) & start(parent.a) <= start(current.a) & end(parent.a) >= end(current.a)) ) { stop("anchor ranges not nested in parent") }
		if (! all(seqnames(parent.t)==seqnames(current.t) & start(parent.t) <= start(current.t) & end(parent.t) >= end(current.t)) ) { stop("target ranges not nested in parent") }
	}

	if (minbox) { 
		# Checking that the minimum bounding box was correctly assigned.
		all.anchors <- all.targets <- list()
		for (x in 1:length(output$indices)) { 
			all.anchors[[x]] <- anchors(collected[[x]])
			all.targets[[x]] <- targets(collected[[x]])
		}
		gathered.a <- unlist(range(split(do.call(c, all.anchors), unlist(output$indices))))
		gathered.t <- unlist(range(split(do.call(c, all.targets), unlist(output$indices))))
		names(gathered.a) <- names(gathered.t) <- NULL
		if (!identical(gathered.a, output$anchors) || !identical(gathered.t, output$targets)) { stop("mismatch in bounding boxes") }
	}	

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
