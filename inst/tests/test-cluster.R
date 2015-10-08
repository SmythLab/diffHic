####################################################################################################
# This tests the clusterPairs function.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

####################################################################################################

simgen <- function(alln, chromos, width, min.space, max.space) {
	# Randomly sampling chromosomes to generate both a pair matrix and the chromosome space.
	output <- GRanges()
	for (x in names(chromos)) {
		n <- chromos[[x]]
		gaps <- round(runif(n, min.space, max.space))
		starts <- cumsum(gaps)
		ends <- starts + width
		suppressWarnings(output <- c(output, GRanges(x, IRanges(starts, ends))))
	}
	
	# Randomly sampling pairs.
   	total <- length(output)
   	chosen1 <- round(runif(alln, 1, total))
   	chosen2 <- round(runif(alln, 1, total))
   	chosen.a <- pmax(chosen1, chosen2)
   	chosen.t <- pmin(chosen1, chosen2)

	# Enforcing uniqueness and anchor/target properties.
   	o <- order(chosen.a, chosen.t)
   	chosen.a <- chosen.a[o]
   	chosen.t <- chosen.t[o]
   	is.diff <- c(TRUE, diff(chosen.a)!=0 | diff(chosen.t)!=0)
	return(DIList(anchors=chosen.a, targets=chosen.t, 
		counts=matrix(0, nrow=alln, ncol=1), totals=0, region=output))
}

crisscross <- function(id1, id2) {
	out <- split(id1, id2)
	if (!all(sapply(out, FUN=function(x) { length(unique(x))==1L }))) {
		return(FALSE)
	}
	out <- split(id2, id1)
	if (!all(sapply(out, FUN=function(x) { length(unique(x))==1L }))) {
		return(FALSE)
	}
	return(TRUE)
}

clustercomp <- function(data, tol, maxw, split=FALSE, data2=NULL) {
	if (!is.null(data2)) { 
		original <- data
		new.regs <- c(regions(data), regions(data2))
		new.as <- c(anchors(data, id=TRUE), anchors(data2, id=TRUE) + length(regions(data)))
		new.ts <- c(targets(data, id=TRUE), targets(data2, id=TRUE) + length(regions(data)))
		o <- order(new.regs)
		reordered <- integer(length(new.regs))
		reordered[o] <- 1:length(new.regs)
		new.as <- reordered[new.as]
		new.ts <- reordered[new.ts]
		new.regs <- new.regs[o]
		data <- DIList(counts=matrix(0, length(new.as), 1), anchors=new.as,
			targets=new.ts, regions=new.regs)
	}

	# Simulating cluster formation first, by expanding each region and checking for overlaps.
	# We use a simple quadratic-time algorithm; slow, but gets the job done.
	region <- regions(data)
	np <- nrow(data)
	expanded <- resize(region, fix="center", width(region)+tol*2)
	allap <- findOverlaps(expanded, region)
#	nested <- findOverlaps(region, region, type="within")
	impossible <- np+1L
	myids <- rep(impossible, np)
	last.id <- 1L

	for (x in 1:np) {
		cura <- data@anchors[x]
		keep.a <- subjectHits(allap)[queryHits(allap)==cura]
		curt <- data@targets[x]
		keep.t <- subjectHits(allap)[queryHits(allap)==curt]

#		has.nested.a<- subjectHits(nested)[queryHits(nested)==cura & subjectHits(nested)!=cura] 
#		has.nested.t <- subjectHits(nested)[queryHits(nested)==curt & subjectHits(nested)!=curt] 
#		anchors.nested <- data@anchors %in% has.nested.a
#		targets.nested <- data@targets %in% has.nested.t
#		if (any(anchors.nested & data@targets %in% keep.t)) { cat("Hooray, nested anchor!\n") }
#		if (any(data@anchors %in% keep.a & targets.nested)) { cat("Hooray, nested target!\n") }
#		if (any(anchors.nested & targets.nested)) { cat("Hooray, nested both!\n") }

		partners <- which(data@anchors %in% keep.a & data@targets %in% keep.t)
		partners <- partners[partners>=x]
		curids <- myids[partners]
		curids <- curids[curids!=impossible]
		
		chosen <- unique(curids)
		if (length(chosen)>1L) { 
			chosen <- min(curids)
			myids[myids%in%curids]<-chosen
		} else if (!length(curids)) {
			chosen <- last.id
			last.id <- last.id + 1L
		} 
		myids[partners]	<- chosen
	}

	if (!is.null(data2)) { 
		comp <- clusterPairs(original, data2, tol=tol, upper=NULL)$indices
		comp <- unlist(comp)
	} else {
		comp <- clusterPairs(data, tol=tol, upper=NULL)$indices[[1]]
	} 
	if (!crisscross(comp, myids)) { stop("mismatches in cluster IDs without bin size restriction") }

# Some error checking functions; just define 'current' as the pairs of interest.
# targets <- which(sapply(split(comp, myids), FUN=function(x) { length(unique(x))!=1 }))
# current <- pairs[myid==targets[1],]
# plot(0,0,xlim=c(1075, 1358),ylim=c(1163,1475))
# rect(start(region)[current$t], start(region)[current$t], end(region)[current$a], end(region)[current$t], col=rgb(1,0,0,0.5))

	# Now, splitting each cluster to keep them under maxw.
	clusters <- split(1:np, comp)
	all.starts <- start(region)
	all.ends <- end(region)+1L
	last <- 0
	myid2 <- myids
	for (x in names(clusters)) {
		active <- clusters[[x]]
		active.a <- data@anchors[active]
		active.t <- data@targets[active]

		cluster.as <- min(all.starts[active.a])
		cluster.ae <- max(all.ends[active.a])
		cluster.ts <- min(all.starts[active.t])
		cluster.te <- max(all.ends[active.t])

		diff.a <- cluster.ae - cluster.as 
		mult.a <- max(1, round(diff.a / maxw))
		jump.a <- diff.a/mult.a

		diff.t <- cluster.te - cluster.ts 
		mult.t <- max(1, round(diff.t / maxw))
		jump.t <- diff.t/mult.t

		mid.a <- (all.starts[active.a]+all.ends[active.a]) * 0.5
		ax <- floor((mid.a - cluster.as)/jump.a)
		mid.t <- (all.starts[active.t]+all.ends[active.t]) * 0.5
		tx <- floor((mid.t - cluster.ts)/jump.t)

		myid2[active] <- ax * mult.t + tx + last
		last <- last + mult.t*mult.a
	}

	# Comparing it to the actual clustering.
	if (!is.null(data2)) { 
		comp2 <- clusterPairs(original, data2, tol=tol, upper=maxw)
		id.comp <- unlist(comp2$indices)
	} else {
		comp2 <- clusterPairs(data, tol=tol, upper=maxw)
		id.comp <- comp2$indices[[1]]
	} 
	if (!crisscross(id.comp, myid2)) { stop("mismatches in cluster IDs when a maximum bin size is applied") }

	# Checking the bounding boxes.
	arange <- range(split(anchors(data), myid2))
	trange <- range(split(targets(data), myid2))
	names(arange) <- names(trange) <- NULL
	if (!identical(comp2$anchors, unlist(arange)) || !identical(comp2$targets, unlist(trange))) {
		stop("mismatches in anchor/target bounding box coordinates")
	}

	if (split) { 
		return(summary(tabulate(id.comp)))
	} else {
		return(summary(tabulate(comp)))
	}
}

####################################################################################################

set.seed(3413094)

chromos <- c(chrA=10, chrB=20, chrC=40)
data <- simgen(100, chromos, 20, 50, 100)
clustercomp(data, tol=70, maxw=200)
clustercomp(data, tol=100, maxw=200)
clustercomp(data, tol=200, maxw=200)

data <- simgen(100, chromos, 10, 50, 100)
clustercomp(data, tol=70, maxw=500)
clustercomp(data, tol=100, maxw=500)
clustercomp(data, tol=200, maxw=500)

data <- simgen(500, chromos, 20, 50, 100)
clustercomp(data, tol=70, maxw=200)
clustercomp(data, tol=100, maxw=200)
clustercomp(data, tol=200, maxw=200)

data <- simgen(500, chromos, 10, 50, 100)
clustercomp(data, tol=70, maxw=500)
clustercomp(data, tol=100, maxw=500)
clustercomp(data, tol=200, maxw=500)

data <- simgen(1000, chromos, 20, 50, 100)
clustercomp(data, tol=70, maxw=200)
clustercomp(data, tol=100, maxw=200)
clustercomp(data, tol=200, maxw=200)

data <- simgen(1000, chromos, 10, 50, 100)
clustercomp(data, tol=70, maxw=500)
clustercomp(data, tol=100, maxw=500)
clustercomp(data, tol=200, maxw=500)

# And again, with flipped settings.	
data <- simgen(100, chromos, 50, 10, 20)
clustercomp(data, tol=0, maxw=200)
clustercomp(data, tol=20, maxw=200)
clustercomp(data, tol=50, maxw=200)

data <- simgen(100, chromos, 100, 10, 20)
clustercomp(data, tol=0, maxw=500)
clustercomp(data, tol=20, maxw=500)
clustercomp(data, tol=50, maxw=500)

data <- simgen(500, chromos, 50, 10, 20)
clustercomp(data, tol=0, maxw=100, split=TRUE)
clustercomp(data, tol=0, maxw=200, split=TRUE)
clustercomp(data, tol=0, maxw=500, split=TRUE)

data <- simgen(500, chromos, 100, 10, 20)
clustercomp(data, tol=0, maxw=100, split=TRUE)
clustercomp(data, tol=0, maxw=200, split=TRUE)
clustercomp(data, tol=0, maxw=500, split=TRUE)

data <- simgen(1000, chromos, 50, 10, 20)
clustercomp(data, tol=0, maxw=100, split=TRUE)
clustercomp(data, tol=0, maxw=200, split=TRUE)
clustercomp(data, tol=0, maxw=500, split=TRUE)

data <- simgen(1000, chromos, 100, 10, 20)
clustercomp(data, tol=0, maxw=100, split=TRUE)
clustercomp(data, tol=0, maxw=200, split=TRUE)
clustercomp(data, tol=0, maxw=500, split=TRUE)

####################################################################################################
# Clustering with multiple entries.

set.seed(34133424)

data <- simgen(100, chromos, 20, 50, 100)
data2 <- simgen(100, chromos, 10, 20, 50)
clustercomp(data, tol=70, maxw=200, data2=data2)
clustercomp(data, tol=100, maxw=200, data2=data2)
clustercomp(data, tol=200, maxw=200, data2=data2)

data <- simgen(100, chromos, 50, 50, 100)
data2 <- simgen(100, chromos, 10, 20, 50)
clustercomp(data, tol=70, maxw=500, data2=data2)
clustercomp(data, tol=100, maxw=500, data2=data2)
clustercomp(data, tol=200, maxw=500, data2=data2)

data <- simgen(500, chromos, 20, 50, 100)
data2 <- simgen(500, chromos, 10, 20, 50)
clustercomp(data, tol=70, maxw=200, data2=data2)
clustercomp(data, tol=100, maxw=200, data2=data2)
clustercomp(data, tol=200, maxw=200, data2=data2)

data <- simgen(500, chromos, 50, 50, 100)
data2 <- simgen(500, chromos, 10, 20, 50)
clustercomp(data, tol=70, maxw=500, data2=data2)
clustercomp(data, tol=100, maxw=500, data2=data2)
clustercomp(data, tol=200, maxw=500, data2=data2)

data <- simgen(1000, chromos, 50, 50, 100)
data2 <- simgen(1000, chromos, 10, 20, 50)
clustercomp(data, tol=70, maxw=200, data2=data2)
clustercomp(data, tol=100, maxw=200, data2=data2)
clustercomp(data, tol=200, maxw=200, data2=data2)

data <- simgen(1000, chromos, 20, 50, 100)
data2 <- simgen(1000, chromos, 10, 20, 50)
clustercomp(data, tol=70, maxw=500, data2=data2)
clustercomp(data, tol=100, maxw=500, data2=data2)
clustercomp(data, tol=200, maxw=500, data2=data2)

####################################################################################################
# Specific clustering tests involving nested things.

require(diffHic)
getLandscape <- function(astarts, aends, tstarts, tends) { 
	astarts <- as.integer(astarts)
	aends <- as.integer(aends)
	tstarts <- as.integer(tstarts)
	tends <- as.integer(tends)
	o <- order(astarts, tstarts)
	.Call(diffHic:::cxx_cluster_2d, astarts[o], tstarts[o], aends[o], tends[o], tol=1L, TRUE)
	invisible(NULL)
}

parent.a <- c(10, 20)
parent.t <- c(10, 20)

# Wholly nested.

nest.a <- c(12, 15)
nest.t <- c(12, 15)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

# Nested target, joined anchor

nest.a <- c(18, 21)
nest.t <- c(12, 15)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

nest.t <- c(10, 12)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

nest.t <- c(15, 20)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

# Nested anchor, joined target.
	
nest.t <- c(18, 21)
nest.a <- c(12, 15)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

nest.a <- c(10, 12)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

nest.a <- c(15, 20)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

# Nested anchor, extended target.

nest.a <- c(12, 15)
nest.t <- c(5, 25)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

# Nested target, extended anchor.

nest.a <- c(5, 25)
nest.t <- c(12, 15)
getLandscape(c(parent.a[1], nest.a[1]),	c(parent.a[2], nest.a[2]), 
	c(parent.t[1], nest.t[1]), c(parent.t[2], nest.t[2]))

# Double nested. 

parent.a2 <- c(10, 20)
parent.t2 <- c(30, 40)

nest.a <- c(12, 15)
nest.t <- c(5, 45)
getLandscape(c(parent.a[1], parent.a2[1], nest.a[1]), c(parent.a[2], parent.a2[2], nest.a[2]), 
	c(parent.t[1], parent.t2[1], nest.t[1]), c(parent.t[2], parent.t2[2], nest.t[2]))

nest.t <- c(15, 35)
getLandscape(c(parent.a[1], parent.a2[1], nest.a[1]), c(parent.a[2], parent.a2[2], nest.a[2]), 
	c(parent.t[1], parent.t2[1], nest.t[1]), c(parent.t[2], parent.t2[2], nest.t[2]))

####################################################################################################
# End.
