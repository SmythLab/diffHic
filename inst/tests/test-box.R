###################################################################################################
# This tests the bin square summarization method in diffHic.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
source("simcounts.R")

chromos<-c(chrA=50, chrB=100, chrC=80)
comp <- function(reference, widths) {
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

	output<- do.call(boxPairs, c(collected, reference=reference))
	stopifnot(length(output$indices)==length(widths))
	for (x in 1:length(output$indices)) { 
		curdex <- output$indices[[x]]
		curlist <- collected[[x]]

		# Checking that each bin pair is truly nested within its reported parent.
		parent.a <- anchors(output$pairs)[curdex]
		parent.t <- targets(output$pairs)[curdex]
		current.a <- anchors(curlist)
		current.t <- targets(curlist)
		
		if (! all(start(parent.a) <= start(current.a) & end(parent.a) >= end(current.a)) ) { stop("anchor ranges not nested in parent") }
		if (! all(start(parent.t) <= start(current.t) & end(parent.t) >= end(current.t)) ) { stop("target ranges not nested in parent") }

		# Subtracting it from the counts.
		if (!identical(counts(output$pairs)[,x], tabulate(output$indices[[x]], nbins=nrow(output$pairs)))) { stop("incidence counts don't match up") }
	}

	return(head(data.frame(anchor.id=output$pairs@anchors, target.id=output$pairs@targets)))
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

####################################################################################################
# End.
