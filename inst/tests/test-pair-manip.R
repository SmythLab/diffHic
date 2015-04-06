####################################################################################################
# This script tests the pair manipulation functions of savePairs and mergePairs.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))
suppressPackageStartupMessages(require(rhdf5))

tmp<-"temp-pairs"
dir.create(tmp)
savecomp<-function(n, nfrags, nchrs) {
	# Simulating the dudes (target<=anchor at all times).
	ai<-as.integer(runif(n, 1, nfrags))
	ti<-as.integer(runif(n, 1, nfrags))
	collected <- data.frame(anchor.id=pmax(ai, ti), target.id=pmin(ai, ti), junk=ai+ti, more.junk=ai-ti)

	# Simulating the fragment IDs.
	blah<-GRanges(sample(paste0("chr", 1:nchrs), nfrags, replace=TRUE), IRanges(1:nfrags, 1:nfrags+10),
		seqinfo=Seqinfo(seqnames=paste0("chr", 1:nchrs)))
	blah<-sort(blah)
	newdir<-file.path(tmp, "output")
	savePairs(collected, newdir, pairParam(fragments=blah))

	# Checking if everything makes sense.
	chrs<-as.character(seqnames(blah))
	indices <- h5ls(newdir)
	indices <- indices[indices$otype=="H5I_DATASET",]
	regot <- list()	
	for (x in 1:nrow(indices)) {
		reread<-h5read(newdir, file.path(indices$group[x], indices$name[x]))
		for (y in 1:ncol(reread)) { attributes(reread[,y]) <- NULL }
		regot[[x]] <- reread
		o <- order(reread$anchor.id, reread$target.id)
		stopifnot(all(diff(o)==1L)) 

		uniq.a<-unique(chrs[reread[,1]])
		uniq.t<-unique(chrs[reread[,2]])
		if (length(uniq.a)!=1L || length(uniq.t)!=1L) { stop("file contains more than one combination") }			
		if (basename(indices$group[x])!=uniq.a || indices$name[x]!=uniq.t) { stop("file contains the incorrect combination") }
	}

	# Checking that the stored result is the same.
	regot <- do.call(rbind, regot)
	regot <- regot[order(regot$anchor.id, regot$target.id, regot$junk, regot$more.junk),]
	original <- collected[order(collected$anchor.id, collected$target.id, collected$junk, collected$more.junk),]
	rownames(original) <- rownames(regot) <- NULL 
	stopifnot(identical(original, regot))
	head(regot)
}

set.seed(23192382)
savecomp(100, 10, 5)
savecomp(100, 10, 15)
savecomp(100, 10, 25)
savecomp(10, 100, 5)
savecomp(10, 100, 15)
savecomp(10, 100, 25)
savecomp(50, 50, 5)
savecomp(50, 50, 15)
savecomp(50, 50, 25)

####################################################################################################
# Finally, chekcing the merging algorithms.

mergecomp<-function(nl, n, nfrags, nchrs) {
	blah<-GRanges(sample(paste0("chr", 1:nchrs), nfrags, replace=TRUE), IRanges(1:nfrags, 1:nfrags+10),
		seqinfo=Seqinfo(seqnames=paste0("chr", 1:nchrs)))
	blah<-sort(blah)
	allfiles<-list()
	allcounts<-list()
	for (x in 1:nl) {
		# Simulating the dudes (target<=anchor at all times).
		ai<-as.integer(runif(n, 1, nfrags))
		ti<-as.integer(runif(n, 1, nfrags))
		collected <- data.frame(anchor.id=pmax(ai, ti), target.id=pmin(ai, ti), junk=ai+ti, more.junk=ai-ti)
		allcounts[[x]]<-collected
		allfiles[[x]]<-  file.path(tmp, paste0("output_", x))
		savePairs(collected, allfiles[[x]], pairParam(fragments=blah))
	}
	
	# Comparing the combined with a more brutal merger.		
	allfiles<-unlist(allfiles)
	allcounts<-do.call(rbind, allcounts)
	rdir<-file.path(tmp, "output_ref")
	savePairs(allcounts, rdir, pairParam(fragments=blah))

	mdir<-file.path(tmp, "output_merged")
	mergePairs(allfiles, mdir)

	# Comparing internal objects.
	combodirs<-c(mdir, rdir)
	out <- diffHic:::.loadIndices(combodirs, seqlevels(blah))
	for (x in names(out)) {
		for (y in names(out[[x]])) {
			current<-out[[x]][[y]]
			stopifnot(all(current))

			alpha <- h5read(mdir, file.path(x, y))
			alpha <- alpha[do.call(order, alpha),]
			bravo <- h5read(rdir, file.path(x, y))
			bravo <- bravo[do.call(order, bravo),]
			rownames(alpha) <- rownames(bravo) <- NULL
			stopifnot(identical(alpha, bravo))
		}
	}

	top.hit <- names(out)[1]
	return(head(alpha))
}

mergecomp(2, 100, 10, 5)
mergecomp(3, 100, 10, 15)
mergecomp(4, 100, 10, 25)
mergecomp(3, 10, 100, 5)
mergecomp(4, 10, 100, 15)
mergecomp(2, 10, 100, 25)
mergecomp(4, 50, 50, 5)
mergecomp(2, 50, 50, 15)
mergecomp(3, 50, 50, 25)

####################################################################################################
# Cleaning up.

unlink(tmp, recursive=TRUE)

####################################################################################################


